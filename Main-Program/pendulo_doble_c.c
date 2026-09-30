#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <SDL2/SDL.h>

#define PI 3.14159265358979323846

#define WIDTH 900
#define HEIGHT 700

#define MAX_TRAIL 4000

/* ============================================================
   SIMULADOR DE PENDULO DOBLE en C

   Estado:
       y[0] = theta1
       y[1] = omega1
       y[2] = theta2
       y[3] = omega2

   Angulos medidos desde la vertical hacia abajo.
   ============================================================ */

typedef struct
{
    double l1;
    double l2;

    double m1;
    double m2;

    double g;

    double b1;
    double b2;

} Parametros;


/* ============================================================
   DERIVADAS DEL SISTEMA
   ============================================================ */

void derivadas(
    const double estado[4],
    double d_estado[4],
    const Parametros *p)
{
    double theta1 = estado[0];
    double omega1 = estado[1];

    double theta2 = estado[2];
    double omega2 = estado[3];

    double delta = theta1 - theta2;

    double c = cos(delta);
    double s = sin(delta);


    /* Matriz de masa */

    double M11 =
        (p->m1 + p->m2) *
        p->l1 * p->l1;

    double M12 =
        p->m2 *
        p->l1 *
        p->l2 *
        c;

    double M21 = M12;

    double M22 =
        p->m2 *
        p->l2 *
        p->l2;


    /* Fuerzas generalizadas */

    double rhs1 =
        -p->m2 *
        p->l1 *
        p->l2 *
        s *
        omega2 * omega2

        -(p->m1 + p->m2) *
        p->g *
        p->l1 *
        sin(theta1)

        -p->b1 * omega1;


    double rhs2 =
        +p->m2 *
        p->l1 *
        p->l2 *
        s *
        omega1 * omega1

        -p->m2 *
        p->g *
        p->l2 *
        sin(theta2)

        -p->b2 * omega2;


    /* --------------------------------------------------------
       Resolver:

       | M11 M12 | | alpha1 | = | rhs1 |
       | M21 M22 | | alpha2 |   | rhs2 |

       Sistema 2x2 resuelto analiticamente.
       -------------------------------------------------------- */

    double determinante =
        M11 * M22 -
        M12 * M21;


    double alpha1 =
        (rhs1 * M22 -
         M12 * rhs2)
        / determinante;


    double alpha2 =
        (M11 * rhs2 -
         rhs1 * M21)
        / determinante;


    d_estado[0] = omega1;
    d_estado[1] = alpha1;

    d_estado[2] = omega2;
    d_estado[3] = alpha2;
}


/* ============================================================
   RUNGE-KUTTA DE CUARTO ORDEN
   ============================================================ */

void paso_rk4(
    double estado[4],
    double dt,
    const Parametros *p)
{
    double k1[4];
    double k2[4];
    double k3[4];
    double k4[4];

    double temp[4];


    /* k1 */

    derivadas(estado, k1, p);


    /* k2 */

    for (int i = 0; i < 4; i++)
        temp[i] =
            estado[i] +
            0.5 * dt * k1[i];

    derivadas(temp, k2, p);


    /* k3 */

    for (int i = 0; i < 4; i++)
        temp[i] =
            estado[i] +
            0.5 * dt * k2[i];

    derivadas(temp, k3, p);


    /* k4 */

    for (int i = 0; i < 4; i++)
        temp[i] =
            estado[i] +
            dt * k3[i];

    derivadas(temp, k4, p);


    /* Actualizacion */

    for (int i = 0; i < 4; i++)
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
   ENERGIA MECANICA
   ============================================================ */

double energia_mecanica(
    const double estado[4],
    const Parametros *p)
{
    double th1 = estado[0];
    double om1 = estado[1];

    double th2 = estado[2];
    double om2 = estado[3];


    double T =
        0.5 *
        (p->m1 + p->m2) *
        p->l1 * p->l1 *
        om1 * om1

        +

        0.5 *
        p->m2 *
        p->l2 * p->l2 *
        om2 * om2

        +

        p->m2 *
        p->l1 *
        p->l2 *
        om1 *
        om2 *
        cos(th1 - th2);


    double V =
        -(p->m1 + p->m2) *
        p->g *
        p->l1 *
        cos(th1)

        -

        p->m2 *
        p->g *
        p->l2 *
        cos(th2);


    return T + V;
}


/* ============================================================
   CONVERSION COORDENADAS FISICAS -> PANTALLA
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
        120 -
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

    double t_max;


    printf(
        "\n--- Simulador de Pendulo Doble en C---\n\n"
    );


    /* --------------------------------------------------------
       Condiciones iniciales
       -------------------------------------------------------- */

    printf(
        "Angulo inicial brazo 1 "
        "(grados, ej: 120): "
    );

    scanf("%lf", &angulo1_grados);


    printf(
        "Angulo inicial brazo 2 "
        "(grados, ej: -10): "
    );

    scanf("%lf", &angulo2_grados);


    printf(
        "Velocidad angular inicial 1 "
        "(grados/s, ej: 0): "
    );

    scanf("%lf", &velocidad1_grados);


    printf(
        "Velocidad angular inicial 2 "
        "(grados/s, ej: 0): "
    );

    scanf("%lf", &velocidad2_grados);


    /* --------------------------------------------------------
       Geometria
       -------------------------------------------------------- */

    printf(
        "Longitud brazo 1 "
        "(m, ej: 1.0): "
    );

    scanf("%lf", &p.l1);


    printf(
        "Longitud brazo 2 "
        "(m, ej: 1.0): "
    );

    scanf("%lf", &p.l2);


    /* --------------------------------------------------------
       Parametros fisicos
       -------------------------------------------------------- */

    printf(
        "Masa 1 "
        "(kg, ej: 1.0): "
    );

    scanf("%lf", &p.m1);


    printf(
        "Masa 2 "
        "(kg, ej: 1.0): "
    );

    scanf("%lf", &p.m2);


    printf(
        "Gravedad "
        "(m/s^2, ej: 9.81): "
    );

    scanf("%lf", &p.g);


    printf(
        "Friccion articulacion 1 "
        "(N*m*s/rad, ej: 0.02): "
    );

    scanf("%lf", &p.b1);


    printf(
        "Friccion articulacion 2 "
        "(N*m*s/rad, ej: 0.02): "
    );

    scanf("%lf", &p.b2);


    printf(
        "Tiempo de simulacion "
        "(s, ej: 20): "
    );

    scanf("%lf", &t_max);


    /* --------------------------------------------------------
       Estado inicial
       -------------------------------------------------------- */

    double estado[4];

    estado[0] =
        angulo1_grados *
        PI / 180.0;

    estado[1] =
        velocidad1_grados *
        PI / 180.0;

    estado[2] =
        angulo2_grados *
        PI / 180.0;

    estado[3] =
        velocidad2_grados *
        PI / 180.0;


    /* Mismo dt que el programa Python */

    const double dt = 0.005;


    /* ========================================================
       SDL
       ======================================================== */

    if (SDL_Init(SDL_INIT_VIDEO) != 0)
    {
        printf(
            "Error iniciando SDL: %s\n",
            SDL_GetError()
        );

        return 1;
    }


    SDL_Window *window =
        SDL_CreateWindow(
            "Pendulo doble",
            SDL_WINDOWPOS_CENTERED,
            SDL_WINDOWPOS_CENTERED,
            WIDTH,
            HEIGHT,
            SDL_WINDOW_SHOWN
        );


    if (window == NULL)
    {
        printf(
            "Error creando ventana: %s\n",
            SDL_GetError()
        );

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
        printf(
            "Error creando renderer: %s\n",
            SDL_GetError()
        );

        SDL_DestroyWindow(window);
        SDL_Quit();

        return 1;
    }


    /* --------------------------------------------------------
       Escala de dibujo
       -------------------------------------------------------- */

    double longitud_total =
        p.l1 + p.l2;

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

    double tiempo = 0.0;

    Uint32 anterior =
        SDL_GetTicks();


    while (
        ejecutando &&
        tiempo <= t_max
    )
    {
        SDL_Event evento;


        while (SDL_PollEvent(&evento))
        {
            if (
                evento.type ==
                SDL_QUIT
            )
            {
                ejecutando = 0;
            }


            if (
                evento.type ==
                SDL_KEYDOWN
            )
            {
                if (
                    evento.key.keysym.sym ==
                    SDLK_ESCAPE
                )
                {
                    ejecutando = 0;
                }
            }
        }


        /*
           Para mantener aproximadamente
           30-35 FPS hacemos varios pasos
           numericos antes de dibujar.
        */

        const double tiempo_frame = 0.03;

        int pasos =
            (int)(
                tiempo_frame /
                dt
            );


        for (
            int k = 0;
            k < pasos;
            k++
        )
        {
            paso_rk4(
                estado,
                dt,
                &p
            );

            tiempo += dt;
        }


        /* ====================================================
           POSICIONES
           ==================================================== */

        double theta1 =
            estado[0];

        double theta2 =
            estado[2];


        double x1 =
            p.l1 *
            sin(theta1);

        double y1 =
            -p.l1 *
            cos(theta1);


        double x2 =
            x1 +
            p.l2 *
            sin(theta2);

        double y2 =
            y1 -
            p.l2 *
            cos(theta2);


        /* Coordenadas pantalla */

        int px0 =
            pantalla_x(
                0.0,
                escala
            );

        int py0 =
            pantalla_y(
                0.0,
                escala
            );


        int px1 =
            pantalla_x(
                x1,
                escala
            );

        int py1 =
            pantalla_y(
                y1,
                escala
            );


        int px2 =
            pantalla_x(
                x2,
                escala
            );

        int py2 =
            pantalla_y(
                y2,
                escala
            );


        /* ====================================================
           GUARDAR TRAYECTORIA
           ==================================================== */

        SDL_Point nuevo;

        nuevo.x = px2;
        nuevo.y = py2;


        if (
            numero_puntos <
            MAX_TRAIL
        )
        {
            trayectoria[
                numero_puntos
            ] = nuevo;

            numero_puntos++;
        }
        else
        {
            for (
                int i = 1;
                i < MAX_TRAIL;
                i++
            )
            {
                trayectoria[
                    i - 1
                ] =
                trayectoria[i];
            }

            trayectoria[
                MAX_TRAIL - 1
            ] = nuevo;
        }


        /* ====================================================
           LIMPIAR PANTALLA
           ==================================================== */

        SDL_SetRenderDrawColor(
            renderer,
            245,
            245,
            245,
            255
        );

        SDL_RenderClear(renderer);


        /* ====================================================
           TRAYECTORIA
           ==================================================== */

        SDL_SetRenderDrawColor(
            renderer,
            170,
            170,
            170,
            255
        );


        if (
            numero_puntos >
            1
        )
        {
            SDL_RenderDrawLines(
                renderer,
                trayectoria,
                numero_puntos
            );
        }


        /* ====================================================
           BRAZOS
           ==================================================== */

        SDL_SetRenderDrawColor(
            renderer,
            30,
            30,
            30,
            255
        );


        SDL_RenderDrawLine(
            renderer,
            px0,
            py0,
            px1,
            py1
        );


        SDL_RenderDrawLine(
            renderer,
            px1,
            py1,
            px2,
            py2
        );


        /* ====================================================
           PIVOTE
           ==================================================== */

        SDL_SetRenderDrawColor(
            renderer,
            0,
            0,
            0,
            255
        );

        dibujar_circulo(
            renderer,
            px0,
            py0,
            6
        );


        /* ====================================================
           MASA 1
           ==================================================== */

        SDL_SetRenderDrawColor(
            renderer,
            60,
            110,
            220,
            255
        );

        dibujar_circulo(
            renderer,
            px1,
            py1,
            11
        );


        /* ====================================================
           MASA 2
           ==================================================== */

        SDL_SetRenderDrawColor(
            renderer,
            220,
            70,
            70,
            255
        );

        dibujar_circulo(
            renderer,
            px2,
            py2,
            11
        );


        SDL_RenderPresent(
            renderer
        );


        /* ====================================================
           DATOS EN TERMINAL
           ==================================================== */

        static int contador = 0;

        contador++;


        if (
            contador >= 15
        )
        {
            contador = 0;


            double energia =
                energia_mecanica(
                    estado,
                    &p
                );


            printf(
                "\n"
                "t = %6.2f s | "
                "theta1 = %8.2f deg | "
                "theta2 = %8.2f deg | "
                "E = %10.4f J",
                tiempo,
                estado[0] *
                180.0 / PI,
                estado[2] *
                180.0 / PI,
                energia
            );


            fflush(stdout);
        }


        /* Limitar aproximadamente
           a 30-35 FPS */

        Uint32 ahora =
            SDL_GetTicks();

        Uint32 transcurrido =
            ahora - anterior;


        if (
            transcurrido <
            30
        )
        {
            SDL_Delay(
                30 -
                transcurrido
            );
        }


        anterior =
            SDL_GetTicks();
    }


    printf(
        "\n\nSimulacion terminada.\n"
    );


    SDL_DestroyRenderer(
        renderer
    );

    SDL_DestroyWindow(
        window
    );

    SDL_Quit();


    return 0;
}
