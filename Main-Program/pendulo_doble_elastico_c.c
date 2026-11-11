#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <SDL2/SDL.h>

#define PI 3.14159265358979323846

#define WIDTH 900
#define HEIGHT 700

#define MAX_TRAIL 4000

/* ============================================================
   SIMULADOR DE PENDULO DOBLE ELASTICO en C

   Cada brazo es un hilo de seda de arana (seda de "dragline" de
   Nephila, valores tipicos aproximados):

       Modulo de Young      E     ~ 10 GPa
       Resistencia a tension sigma ~ 1.1 GPa

   Con la seccion A de cada brazo:

       rigidez       k  = E * A / L0
       tension max.  Tmax = sigma * A

   Si la tension de un brazo supera Tmax el brazo se rompe, el
   programa termina y avisa cual pendulo cayo.

   Como el brazo ya no es rigido, el estado se escribe en
   coordenadas cartesianas (el pivote esta en el origen, y hacia
   arriba):

       estado[0] = x1    estado[1] = y1
       estado[2] = x2    estado[3] = y2
       estado[4] = vx1   estado[5] = vy1
       estado[6] = vx2   estado[7] = vy2

   Los angulos de entrada se miden desde la vertical hacia abajo,
   igual que antes.
   ============================================================ */

/* Propiedades de la seda */

#define E_SEDA        10.0e9    /* Pa  (modulo de Young)            */
#define SIGMA_SEDA     1.1e9    /* Pa  (resistencia a la tension)   */

/* Amortiguamiento interno del hilo, como fraccion del critico.
   La seda real disipa mucha energia; esto tambien apaga las
   vibraciones axiales rapidas. */

#define ZETA_AXIAL     0.05

/* 1 = el hilo solo tira (se afloja si se comprime, como una
   cuerda).  0 = se comporta como un resorte que tambien empuja. */

#ifndef SOLO_TENSION
#define SOLO_TENSION   1
#endif


typedef struct
{
    double l1;
    double l2;

    double m1;
    double m2;

    double g;

    double b1;
    double b2;

    double area;        /* m^2 */

    double k1;          /* N/m */
    double k2;

    double c1;          /* N*s/m (amortiguamiento axial) */
    double c2;

    double T_rotura;    /* N, igual para ambos brazos */

} Parametros;


/* ============================================================
   FUERZAS Y TENSIONES

   Calcula la tension de cada brazo y la fuerza total sobre
   cada masa.
   ============================================================ */

void calcular(
    const double e[8],
    const Parametros *p,
    double *T1,
    double *T2,
    double F[4])
{
    double x1 = e[0];
    double y1 = e[1];
    double x2 = e[2];
    double y2 = e[3];

    double vx1 = e[4];
    double vy1 = e[5];
    double vx2 = e[6];
    double vy2 = e[7];


    /* Brazo 1: pivote -> masa 1 */

    double dx1 = x1;
    double dy1 = y1;

    double L1 = sqrt(dx1 * dx1 + dy1 * dy1);
    if (L1 < 1e-9) L1 = 1e-9;

    double u1x = dx1 / L1;
    double u1y = dy1 / L1;


    /* Brazo 2: masa 1 -> masa 2 */

    double dx2 = x2 - x1;
    double dy2 = y2 - y1;

    double L2 = sqrt(dx2 * dx2 + dy2 * dy2);
    if (L2 < 1e-9) L2 = 1e-9;

    double u2x = dx2 / L2;
    double u2y = dy2 / L2;


    /* Velocidad de alargamiento de cada brazo */

    double vr1 = vx1 * u1x + vy1 * u1y;

    double vr2 =
        (vx2 - vx1) * u2x +
        (vy2 - vy1) * u2y;


    /* Tension (positiva = el hilo tira) */

    double ext1 = L1 - p->l1;
    double ext2 = L2 - p->l2;

    double t1 = p->k1 * ext1 + p->c1 * vr1;
    double t2 = p->k2 * ext2 + p->c2 * vr2;

#if SOLO_TENSION
    if (ext1 <= 0.0 || t1 < 0.0) t1 = 0.0;
    if (ext2 <= 0.0 || t2 < 0.0) t2 = 0.0;
#endif

    *T1 = t1;
    *T2 = t2;


    /* Fuerzas: gravedad + tension */

    double fx1 = -t1 * u1x + t2 * u2x;
    double fy1 = -p->m1 * p->g - t1 * u1y + t2 * u2y;

    double fx2 = -t2 * u2x;
    double fy2 = -p->m2 * p->g - t2 * u2y;


    /* Friccion en las articulaciones: torque -b * omega, con
       omega la velocidad angular absoluta de cada brazo.
       Se aplica como fuerza tangencial sobre la masa del brazo. */

    double w1 = (x1 * vy1 - y1 * vx1) / (L1 * L1);

    double w2 =
        (dx2 * (vy2 - vy1) - dy2 * (vx2 - vx1)) /
        (L2 * L2);

    double ft1 = -p->b1 * w1 / L1;
    double ft2 = -p->b2 * w2 / L2;

    fx1 += ft1 * (-dy1 / L1);
    fy1 += ft1 * ( dx1 / L1);

    fx2 += ft2 * (-dy2 / L2);
    fy2 += ft2 * ( dx2 / L2);


    F[0] = fx1;
    F[1] = fy1;
    F[2] = fx2;
    F[3] = fy2;
}


/* ============================================================
   DERIVADAS DEL SISTEMA
   ============================================================ */

void derivadas(
    const double estado[8],
    double d_estado[8],
    const Parametros *p)
{
    double T1, T2;
    double F[4];

    calcular(estado, p, &T1, &T2, F);

    d_estado[0] = estado[4];
    d_estado[1] = estado[5];
    d_estado[2] = estado[6];
    d_estado[3] = estado[7];

    d_estado[4] = F[0] / p->m1;
    d_estado[5] = F[1] / p->m1;
    d_estado[6] = F[2] / p->m2;
    d_estado[7] = F[3] / p->m2;
}


/* ============================================================
   RUNGE-KUTTA DE CUARTO ORDEN
   ============================================================ */

void paso_rk4(
    double estado[8],
    double dt,
    const Parametros *p)
{
    double k1[8];
    double k2[8];
    double k3[8];
    double k4[8];

    double temp[8];


    derivadas(estado, k1, p);

    for (int i = 0; i < 8; i++)
        temp[i] = estado[i] + 0.5 * dt * k1[i];

    derivadas(temp, k2, p);

    for (int i = 0; i < 8; i++)
        temp[i] = estado[i] + 0.5 * dt * k2[i];

    derivadas(temp, k3, p);

    for (int i = 0; i < 8; i++)
        temp[i] = estado[i] + dt * k3[i];

    derivadas(temp, k4, p);


    for (int i = 0; i < 8; i++)
    {
        estado[i] +=
            (dt / 6.0) *
            (
                k1[i]
                + 2.0 * k2[i]
                + 2.0 * k3[i]
                + k4[i]
            );
    }
}


/* ============================================================
   ENERGIA MECANICA (cinetica + gravitatoria + elastica)
   ============================================================ */

double energia_mecanica(
    const double e[8],
    const Parametros *p)
{
    double x1 = e[0];
    double y1 = e[1];
    double x2 = e[2];
    double y2 = e[3];

    double T =
        0.5 * p->m1 * (e[4] * e[4] + e[5] * e[5]) +
        0.5 * p->m2 * (e[6] * e[6] + e[7] * e[7]);

    double V =
        p->m1 * p->g * y1 +
        p->m2 * p->g * y2;

    double L1 = sqrt(x1 * x1 + y1 * y1);
    double L2 = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));

    double ext1 = L1 - p->l1;
    double ext2 = L2 - p->l2;

#if SOLO_TENSION
    if (ext1 < 0.0) ext1 = 0.0;
    if (ext2 < 0.0) ext2 = 0.0;
#endif

    V += 0.5 * p->k1 * ext1 * ext1;
    V += 0.5 * p->k2 * ext2 * ext2;

    return T + V;
}


/* ============================================================
   CONVERSION COORDENADAS FISICAS -> PANTALLA
   (pivote en el centro de la ventana)
   ============================================================ */

int pantalla_x(
    double x,
    double escala)
{
    return
        WIDTH / 2 +
        (int)(x * escala);
}


int pantalla_y(
    double y,
    double escala)
{
    return
        HEIGHT / 2 -
        (int)(y * escala);
}


/* ============================================================
   CIRCULO SOLIDO
   ============================================================ */

void dibujar_circulo(
    SDL_Renderer *renderer,
    int cx,
    int cy,
    int radio)
{
    for (int y = -radio; y <= radio; y++)
    {
        for (int x = -radio; x <= radio; x++)
        {
            if (x*x + y*y <= radio*radio)
            {
                SDL_RenderDrawPoint(
                    renderer,
                    cx + x,
                    cy + y
                );
            }
        }
    }
}


/* ============================================================
   MAIN
   ============================================================ */

int main(void)
{
    Parametros p;

    double angulo1_grados;
    double angulo2_grados;

    double velocidad1_grados;
    double velocidad2_grados;

    double area_cm2;

    double t_max;


    printf(
        "\n--- Simulador de Pendulo Doble Elastico (seda de arana) en C ---\n\n"
    );


    /* --------------------------------------------------------
       Condiciones iniciales
       -------------------------------------------------------- */

    printf("Angulo inicial brazo 1 (grados, ej: 120): ");
    scanf("%lf", &angulo1_grados);

    printf("Angulo inicial brazo 2 (grados, ej: -10): ");
    scanf("%lf", &angulo2_grados);

    printf("Velocidad angular inicial 1 (grados/s, ej: 0): ");
    scanf("%lf", &velocidad1_grados);

    printf("Velocidad angular inicial 2 (grados/s, ej: 0): ");
    scanf("%lf", &velocidad2_grados);


    /* --------------------------------------------------------
       Geometria
       -------------------------------------------------------- */

    printf("Longitud brazo 1 (m, ej: 1.0): ");
    scanf("%lf", &p.l1);

    printf("Longitud brazo 2 (m, ej: 1.0): ");
    scanf("%lf", &p.l2);

    printf("Seccion de cada brazo (cm^2, ej: 15): ");
    scanf("%lf", &area_cm2);


    /* --------------------------------------------------------
       Parametros fisicos
       -------------------------------------------------------- */

    printf("Masa 1 (kg, ej: 1.0): ");
    scanf("%lf", &p.m1);

    printf("Masa 2 (kg, ej: 1.0): ");
    scanf("%lf", &p.m2);

    printf("Gravedad (m/s^2, ej: 9.81): ");
    scanf("%lf", &p.g);

    printf("Friccion articulacion 1 (N*m*s/rad, ej: 0.02): ");
    scanf("%lf", &p.b1);

    printf("Friccion articulacion 2 (N*m*s/rad, ej: 0.02): ");
    scanf("%lf", &p.b2);

    printf("Tiempo de simulacion (s, ej: 20): ");
    scanf("%lf", &t_max);


    if (p.l1 <= 0.0 || p.l2 <= 0.0 || p.m1 <= 0.0 || p.m2 <= 0.0 ||
        area_cm2 <= 0.0)
    {
        printf("Longitudes, masas y seccion deben ser positivas.\n");
        return 1;
    }


    /* --------------------------------------------------------
       Propiedades derivadas de la seda
       -------------------------------------------------------- */

    p.area = area_cm2 * 1.0e-4;                 /* cm^2 -> m^2 */

    p.k1 = E_SEDA * p.area / p.l1;
    p.k2 = E_SEDA * p.area / p.l2;

    p.c1 = 2.0 * ZETA_AXIAL * sqrt(p.k1 * (p.m1 + p.m2));
    p.c2 = 2.0 * ZETA_AXIAL * sqrt(p.k2 * p.m2);

    p.T_rotura = SIGMA_SEDA * p.area;

    printf(
        "\nSeda: E = %.1f GPa, resistencia = %.2f GPa\n"
        "Rigidez  k1 = %.3e N/m, k2 = %.3e N/m\n"
        "Cada brazo se rompe con una tension de %.4g N (~%.4g kgf)\n",
        E_SEDA / 1e9,
        SIGMA_SEDA / 1e9,
        p.k1,
        p.k2,
        p.T_rotura,
        p.T_rotura / 9.81
    );


    /* --------------------------------------------------------
       Estado inicial (brazos sin estirar)
       -------------------------------------------------------- */

    double th1 = angulo1_grados * PI / 180.0;
    double th2 = angulo2_grados * PI / 180.0;

    double om1 = velocidad1_grados * PI / 180.0;
    double om2 = velocidad2_grados * PI / 180.0;

    double estado[8];

    estado[0] =  p.l1 * sin(th1);
    estado[1] = -p.l1 * cos(th1);

    estado[2] = estado[0] + p.l2 * sin(th2);
    estado[3] = estado[1] - p.l2 * cos(th2);

    estado[4] = om1 * p.l1 * cos(th1);
    estado[5] = om1 * p.l1 * sin(th1);

    estado[6] = estado[4] + om2 * p.l2 * cos(th2);
    estado[7] = estado[5] + om2 * p.l2 * sin(th2);


    /* --------------------------------------------------------
       Paso de integracion.

       La seda es muy rigida: su frecuencia propia es de miles de
       rad/s, asi que el paso debe ser mucho menor que el del
       pendulo rigido (~ 0.1 / omega_max, y nunca mas de 0.005).
       -------------------------------------------------------- */

    double omega2_max =
        (p.k1 + 2.0 * p.k2) * (1.0 / p.m1 + 1.0 / p.m2);

    double dt_max = 0.1 / sqrt(omega2_max);

    if (dt_max > 0.005)
        dt_max = 0.005;

    const double tiempo_frame = 0.03;

    int pasos = (int)ceil(tiempo_frame / dt_max);

    const double dt = tiempo_frame / pasos;

    printf(
        "Paso de integracion dt = %.3e s (%d pasos RK4 por cuadro)\n",
        dt,
        pasos
    );

    if (pasos > 20000)
        printf(
            "Aviso: con esta rigidez la simulacion correra mas lenta "
            "que el tiempo real.\n"
        );


    /* ========================================================
       SDL
       ======================================================== */

    if (SDL_Init(SDL_INIT_VIDEO) != 0)
    {
        printf("Error iniciando SDL: %s\n", SDL_GetError());
        return 1;
    }


    SDL_Window *window =
        SDL_CreateWindow(
            "Pendulo doble elastico",
            SDL_WINDOWPOS_CENTERED,
            SDL_WINDOWPOS_CENTERED,
            WIDTH,
            HEIGHT,
            SDL_WINDOW_SHOWN
        );

    if (window == NULL)
    {
        printf("Error creando ventana: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }


    SDL_Renderer *renderer =
        SDL_CreateRenderer(
            window,
            -1,
            SDL_RENDERER_ACCELERATED |
            SDL_RENDERER_PRESENTVSYNC
        );

    if (renderer == NULL)
    {
        printf("Error creando renderer: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }


    /* --------------------------------------------------------
       Escala de dibujo
       -------------------------------------------------------- */

    double longitud_total = p.l1 + p.l2;

    double escala =
        0.75 *
        (HEIGHT / 2.0) /
        longitud_total;


    /* --------------------------------------------------------
       Trayectoria masa 2
       -------------------------------------------------------- */

    SDL_Point trayectoria[MAX_TRAIL];

    int numero_puntos = 0;


    /* ========================================================
       BUCLE PRINCIPAL
       ======================================================== */

    int ejecutando = 1;

    int roto1 = 0;
    int roto2 = 0;

    double T1 = 0.0;
    double T2 = 0.0;

    double tiempo = 0.0;

    Uint32 anterior = SDL_GetTicks();


    while (
        ejecutando &&
        tiempo <= t_max
    )
    {
        SDL_Event evento;

        while (SDL_PollEvent(&evento))
        {
            if (evento.type == SDL_QUIT)
                ejecutando = 0;

            if (evento.type == SDL_KEYDOWN)
            {
                if (evento.key.keysym.sym == SDLK_ESCAPE)
                    ejecutando = 0;
            }
        }


        /* ----------------------------------------------------
           Pasos numericos (se revisa la tension en CADA paso,
           porque el tiron del hilo al tensarse dura microsegundos)
           ---------------------------------------------------- */

        for (int k = 0; k < pasos; k++)
        {
            paso_rk4(estado, dt, &p);

            tiempo += dt;

            double F[4];

            calcular(estado, &p, &T1, &T2, F);

            roto1 = (T1 > p.T_rotura);
            roto2 = (T2 > p.T_rotura);

            if (roto1 || roto2)
                break;

            if (!isfinite(estado[0]) || !isfinite(estado[2]))
            {
                printf("\n\nLa simulacion divergio numericamente.\n");
                ejecutando = 0;
                break;
            }
        }


        /* ====================================================
           POSICIONES
           ==================================================== */

        double x1 = estado[0];
        double y1 = estado[1];

        double x2 = estado[2];
        double y2 = estado[3];


        /* Coordenadas pantalla */

        int px0 = pantalla_x(0.0, escala);
        int py0 = pantalla_y(0.0, escala);

        int px1 = pantalla_x(x1, escala);
        int py1 = pantalla_y(y1, escala);

        int px2 = pantalla_x(x2, escala);
        int py2 = pantalla_y(y2, escala);


        /* ====================================================
           GUARDAR TRAYECTORIA
           ==================================================== */

        SDL_Point nuevo;

        nuevo.x = px2;
        nuevo.y = py2;


        if (numero_puntos < MAX_TRAIL)
        {
            trayectoria[numero_puntos] = nuevo;
            numero_puntos++;
        }
        else
        {
            for (int i = 1; i < MAX_TRAIL; i++)
                trayectoria[i - 1] = trayectoria[i];

            trayectoria[MAX_TRAIL - 1] = nuevo;
        }


        /* ====================================================
           LIMPIAR PANTALLA
           ==================================================== */

        SDL_SetRenderDrawColor(renderer, 245, 245, 245, 255);
        SDL_RenderClear(renderer);


        /* ====================================================
           TRAYECTORIA
           ==================================================== */

        SDL_SetRenderDrawColor(renderer, 170, 170, 170, 255);

        if (numero_puntos > 1)
            SDL_RenderDrawLines(renderer, trayectoria, numero_puntos);


        /* ====================================================
           BRAZOS (el que se rompio se dibuja en rojo)
           ==================================================== */

        if (roto1)
            SDL_SetRenderDrawColor(renderer, 220, 40, 40, 255);
        else
            SDL_SetRenderDrawColor(renderer, 30, 30, 30, 255);

        SDL_RenderDrawLine(renderer, px0, py0, px1, py1);


        if (roto2)
            SDL_SetRenderDrawColor(renderer, 220, 40, 40, 255);
        else
            SDL_SetRenderDrawColor(renderer, 30, 30, 30, 255);

        SDL_RenderDrawLine(renderer, px1, py1, px2, py2);


        /* ====================================================
           PIVOTE
           ==================================================== */

        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        dibujar_circulo(renderer, px0, py0, 6);


        /* ====================================================
           MASA 1
           ==================================================== */

        SDL_SetRenderDrawColor(renderer, 60, 110, 220, 255);
        dibujar_circulo(renderer, px1, py1, 11);


        /* ====================================================
           MASA 2
           ==================================================== */

        SDL_SetRenderDrawColor(renderer, 220, 70, 70, 255);
        dibujar_circulo(renderer, px2, py2, 11);


        SDL_RenderPresent(renderer);


        /* ====================================================
           ROTURA: avisar y terminar
           ==================================================== */

        if (roto1 || roto2)
        {
            char mensaje[512];

            if (roto1 && roto2)
                snprintf(
                    mensaje, sizeof(mensaje),
                    "Se rompieron los brazos 1 y 2 a t = %.4f s.\n"
                    "Cayeron el PENDULO 1 y el PENDULO 2.\n"
                    "Tension brazo 1: %.4g N, brazo 2: %.4g N "
                    "(limite %.4g N).",
                    tiempo, T1, T2, p.T_rotura
                );
            else if (roto1)
                snprintf(
                    mensaje, sizeof(mensaje),
                    "Se rompio el brazo 1 a t = %.4f s.\n"
                    "Cayo el PENDULO 1 (y con el, el pendulo 2).\n"
                    "Tension: %.4g N (limite %.4g N).",
                    tiempo, T1, p.T_rotura
                );
            else
                snprintf(
                    mensaje, sizeof(mensaje),
                    "Se rompio el brazo 2 a t = %.4f s.\n"
                    "Cayo el PENDULO 2.\n"
                    "Tension: %.4g N (limite %.4g N).",
                    tiempo, T2, p.T_rotura
                );

            printf("\n\n*** %s ***\n", mensaje);
            fflush(stdout);

            SDL_ShowSimpleMessageBox(
                SDL_MESSAGEBOX_WARNING,
                "Brazo roto",
                mensaje,
                window
            );

            break;
        }


        /* ====================================================
           DATOS EN TERMINAL
           ==================================================== */

        static int contador = 0;

        contador++;

        if (contador >= 15)
        {
            contador = 0;

            double energia = energia_mecanica(estado, &p);

            double th1_deg =
                atan2(estado[0], -estado[1]) * 180.0 / PI;

            double th2_deg =
                atan2(estado[2] - estado[0], -(estado[3] - estado[1])) *
                180.0 / PI;

            printf(
                "\nt = %6.2f s | th1 = %8.2f | th2 = %8.2f | "
                "T1 = %10.4g N (%5.1f%%) | T2 = %10.4g N (%5.1f%%) | "
                "E = %10.4f J",
                tiempo,
                th1_deg,
                th2_deg,
                T1, 100.0 * T1 / p.T_rotura,
                T2, 100.0 * T2 / p.T_rotura,
                energia
            );

            fflush(stdout);
        }


        /* Limitar aproximadamente a 30-35 FPS */

        Uint32 ahora = SDL_GetTicks();
        Uint32 transcurrido = ahora - anterior;

        if (transcurrido < 30)
            SDL_Delay(30 - transcurrido);

        anterior = SDL_GetTicks();
    }


    if (!(roto1 || roto2))
        printf("\n\nSimulacion terminada. Ningun brazo se rompio.\n");


    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();


    return 0;
}
