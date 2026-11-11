#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <SDL2/SDL.h>

#define PI 3.14159265358979323846

#define WIDTH 900
#define HEIGHT 700

#define MAX_TRAIL 4000
#define MAX_N 300

/* ============================================================
   SIMULADOR DE PENDULO N-PLE en C

   Cadena de N masas puntuales unidas por barras rigidas.

   Estado (tamano 2N):
       estado[0 .. N-1]   = theta_i  (angulo desde la vertical hacia abajo)
       estado[N .. 2N-1]  = omega_i

   Ecuacion de movimiento (Lagrange):

       M(theta) * alpha = rhs(theta, omega)

       M_ij   = S_max(i,j) * l_i * l_j * cos(theta_i - theta_j)
       rhs_i  = - sum_j S_max(i,j) l_i l_j sin(theta_i - theta_j) omega_j^2
                - S_i * g * l_i * sin(theta_i)
                - b_i * omega_i

   donde S_i = m_i + m_(i+1) + ... + m_N  (masa "colgando" desde i).

   M es simetrica definida positiva -> se resuelve con Cholesky.
   Con N=2 reproduce exactamente el pendulo doble anterior.
   ============================================================ */

typedef struct
{
    int n;

    double g;

    double *l;      /* longitudes            */
    double *m;      /* masas                 */
    double *S;      /* sumas de masa de i..N */
    double *b;      /* friccion articulacion */

    /* Memoria de trabajo (evita malloc en cada paso) */
    double *M;      /* matriz N x N          */
    double *rhs;    /* vector N              */
    double *y;      /* vector auxiliar N     */

    double *k1, *k2, *k3, *k4, *temp;   /* RK4, tamano 2N */

} Sistema;


int sistema_crear(
    Sistema *s,
    int n,
    double longitud_total,
    double masa,
    double g,
    double friccion)
{
    s->n = n;
    s->g = g;

    s->l    = malloc(n * sizeof(double));
    s->m    = malloc(n * sizeof(double));
    s->S    = malloc(n * sizeof(double));
    s->b    = malloc(n * sizeof(double));

    s->M    = malloc((size_t)n * n * sizeof(double));
    s->rhs  = malloc(n * sizeof(double));
    s->y    = malloc(n * sizeof(double));

    s->k1   = malloc(2 * n * sizeof(double));
    s->k2   = malloc(2 * n * sizeof(double));
    s->k3   = malloc(2 * n * sizeof(double));
    s->k4   = malloc(2 * n * sizeof(double));
    s->temp = malloc(2 * n * sizeof(double));

    if (!s->l || !s->m || !s->S || !s->b || !s->M || !s->rhs ||
        !s->y || !s->k1 || !s->k2 || !s->k3 || !s->k4 || !s->temp)
        return 0;

    for (int i = 0; i < n; i++)
    {
        s->l[i] = longitud_total / n;
        s->m[i] = masa;
        s->b[i] = friccion;
    }

    double acumulada = 0.0;

    for (int i = n - 1; i >= 0; i--)
    {
        acumulada += s->m[i];
        s->S[i] = acumulada;
    }

    return 1;
}


void sistema_liberar(Sistema *s)
{
    free(s->l);  free(s->m);  free(s->S);  free(s->b);
    free(s->M);  free(s->rhs); free(s->y);
    free(s->k1); free(s->k2); free(s->k3); free(s->k4); free(s->temp);
}


/* ============================================================
   DERIVADAS DEL SISTEMA
   ============================================================ */

void derivadas(
    const double *estado,
    double *d_estado,
    Sistema *s)
{
    const int n = s->n;

    const double *th = estado;
    const double *om = estado + n;

    double *M   = s->M;
    double *rhs = s->rhs;


    /* Parte gravitatoria y friccion de rhs */

    for (int i = 0; i < n; i++)
    {
        rhs[i] =
            - s->S[i] * s->g * s->l[i] * sin(th[i])
            - s->b[i] * om[i];
    }


    /* Matriz de masa y terminos centrifugos.
       Recorremos solo j <= i y usamos simetria. */

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            double A = s->S[i] * s->l[i] * s->l[j];   /* max(i,j) = i */

            double delta = th[i] - th[j];
            double c = cos(delta);

            M[i * n + j] = A * c;
            M[j * n + i] = A * c;

            if (j != i)
            {
                double sd = sin(delta);

                rhs[i] -= A * sd * om[j] * om[j];
                rhs[j] += A * sd * om[i] * om[i];
            }
        }
    }


    /* --------------------------------------------------------
       Cholesky in-place (triangulo inferior): M = L * L^T
       -------------------------------------------------------- */

    for (int j = 0; j < n; j++)
    {
        double suma = M[j * n + j];

        for (int k = 0; k < j; k++)
            suma -= M[j * n + k] * M[j * n + k];

        double Ljj = sqrt(suma);
        M[j * n + j] = Ljj;

        for (int i = j + 1; i < n; i++)
        {
            double v = M[i * n + j];

            for (int k = 0; k < j; k++)
                v -= M[i * n + k] * M[j * n + k];

            M[i * n + j] = v / Ljj;
        }
    }


    /* Sustitucion hacia adelante: L y = rhs */

    for (int i = 0; i < n; i++)
    {
        double v = rhs[i];

        for (int k = 0; k < i; k++)
            v -= M[i * n + k] * s->y[k];

        s->y[i] = v / M[i * n + i];
    }


    /* Sustitucion hacia atras: L^T alpha = y */

    for (int i = n - 1; i >= 0; i--)
    {
        double v = s->y[i];

        for (int k = i + 1; k < n; k++)
            v -= M[k * n + i] * d_estado[n + k];

        d_estado[n + i] = v / M[i * n + i];
    }


    for (int i = 0; i < n; i++)
        d_estado[i] = om[i];
}


/* ============================================================
   RUNGE-KUTTA DE CUARTO ORDEN
   ============================================================ */

void paso_rk4(
    double *estado,
    double dt,
    Sistema *s)
{
    const int N2 = 2 * s->n;

    derivadas(estado, s->k1, s);

    for (int i = 0; i < N2; i++)
        s->temp[i] = estado[i] + 0.5 * dt * s->k1[i];

    derivadas(s->temp, s->k2, s);

    for (int i = 0; i < N2; i++)
        s->temp[i] = estado[i] + 0.5 * dt * s->k2[i];

    derivadas(s->temp, s->k3, s);

    for (int i = 0; i < N2; i++)
        s->temp[i] = estado[i] + dt * s->k3[i];

    derivadas(s->temp, s->k4, s);

    for (int i = 0; i < N2; i++)
    {
        estado[i] +=
            (dt / 6.0) *
            (
                s->k1[i]
                + 2.0 * s->k2[i]
                + 2.0 * s->k3[i]
                + s->k4[i]
            );
    }
}


/* ============================================================
   ENERGIA MECANICA:  T = 1/2 w^T M w ,  V = -g sum S_i l_i cos(theta_i)
   ============================================================ */

double energia_mecanica(
    const double *estado,
    const Sistema *s)
{
    const int n = s->n;

    const double *th = estado;
    const double *om = estado + n;

    double T = 0.0;
    double V = 0.0;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            int mx = (i > j) ? i : j;

            T += 0.5 *
                 s->S[mx] * s->l[i] * s->l[j] *
                 cos(th[i] - th[j]) *
                 om[i] * om[j];
        }

        V -= s->S[i] * s->g * s->l[i] * cos(th[i]);
    }

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
    int n;

    double angulo1_grados;
    double angulo_resto_grados;
    double velocidad_grados;

    double longitud_total;
    double masa;
    double g;
    double friccion;
    double t_max;


    printf(
        "\n--- Simulador de Pendulo N-ple en C ---\n\n"
    );


    printf("Numero de pendulos N (1-%d, ej: 2): ", MAX_N);
    if (scanf("%d", &n) != 1 || n < 1 || n > MAX_N)
    {
        printf("N invalido.\n");
        return 1;
    }

    printf("Angulo inicial brazo 1 (grados, ej: 120): ");
    scanf("%lf", &angulo1_grados);

    printf("Angulo inicial de los demas brazos (grados, ej: -10): ");
    scanf("%lf", &angulo_resto_grados);

    printf("Velocidad angular inicial de todos (grados/s, ej: 0): ");
    scanf("%lf", &velocidad_grados);

    printf("Longitud TOTAL de la cadena (m, ej: 2.0): ");
    scanf("%lf", &longitud_total);

    printf("Masa de cada pendulo (kg, ej: 1.0): ");
    scanf("%lf", &masa);

    printf("Gravedad (m/s^2, ej: 9.81): ");
    scanf("%lf", &g);

    printf("Friccion en cada articulacion (N*m*s/rad, ej: 0.02): ");
    scanf("%lf", &friccion);

    printf("Tiempo de simulacion (s, ej: 20): ");
    scanf("%lf", &t_max);


    Sistema sistema;

    if (!sistema_crear(&sistema, n, longitud_total, masa, g, friccion))
    {
        printf("No hay memoria suficiente.\n");
        return 1;
    }


    /* Estado inicial */

    double *estado = calloc(2 * n, sizeof(double));
    int *px = malloc((n + 1) * sizeof(int));
    int *py = malloc((n + 1) * sizeof(int));

    if (!estado || !px || !py)
    {
        printf("No hay memoria suficiente.\n");
        return 1;
    }

    for (int i = 0; i < n; i++)
    {
        double ang = (i == 0) ? angulo1_grados : angulo_resto_grados;

        estado[i]     = ang * PI / 180.0;
        estado[n + i] = velocidad_grados * PI / 180.0;
    }


    /* El paso debe encogerse al crecer N (el sistema se vuelve rigido):
       sin esto la energia se dispara y el solver termina en NaN
       (medido: dt=0.005 ya explota con N=50). */

    const double tiempo_frame = 0.03;

    double dt_max = 0.05 / n;
    if (dt_max > 0.005) dt_max = 0.005;

    int pasos = (int)ceil(tiempo_frame / dt_max);
    const double dt = tiempo_frame / pasos;

    if (n > 100)
        printf(
            "\nAviso: con N=%d la simulacion correra mas lenta que "
            "el tiempo real (%d pasos RK4 por cuadro).\n",
            n, pasos
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
            "Pendulo N-ple",
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


    /* Escala: la cadena extendida ocupa 75% del radio de la ventana */

    double escala =
        0.75 *
        (HEIGHT / 2.0) /
        longitud_total;


    /* Tamano de las masas segun N (para que no se amontonen) */

    int radio_masa = 12 - n / 4;
    if (radio_masa < 2)  radio_masa = 2;
    if (radio_masa > 11) radio_masa = 11;


    /* Trayectoria de la ultima masa */

    SDL_Point *trayectoria = malloc(MAX_TRAIL * sizeof(SDL_Point));
    int numero_puntos = 0;


    /* ========================================================
       BUCLE PRINCIPAL
       ======================================================== */

    int ejecutando = 1;
    double tiempo = 0.0;
    Uint32 anterior = SDL_GetTicks();
    int contador = 0;


    while (ejecutando && tiempo <= t_max)
    {
        SDL_Event evento;

        while (SDL_PollEvent(&evento))
        {
            if (evento.type == SDL_QUIT)
                ejecutando = 0;

            if (evento.type == SDL_KEYDOWN &&
                evento.key.keysym.sym == SDLK_ESCAPE)
                ejecutando = 0;
        }


        for (int k = 0; k < pasos; k++)
        {
            paso_rk4(estado, dt, &sistema);
            tiempo += dt;
        }


        /* ====================================================
           POSICIONES
           ==================================================== */

        double x = 0.0;
        double y = 0.0;

        px[0] = pantalla_x(0.0, escala);
        py[0] = pantalla_y(0.0, escala);

        for (int i = 0; i < n; i++)
        {
            x += sistema.l[i] * sin(estado[i]);
            y -= sistema.l[i] * cos(estado[i]);

            px[i + 1] = pantalla_x(x, escala);
            py[i + 1] = pantalla_y(y, escala);
        }


        /* ====================================================
           GUARDAR TRAYECTORIA (ultima masa)
           ==================================================== */

        SDL_Point nuevo;
        nuevo.x = px[n];
        nuevo.y = py[n];

        if (numero_puntos < MAX_TRAIL)
        {
            trayectoria[numero_puntos++] = nuevo;
        }
        else
        {
            for (int i = 1; i < MAX_TRAIL; i++)
                trayectoria[i - 1] = trayectoria[i];

            trayectoria[MAX_TRAIL - 1] = nuevo;
        }


        /* ====================================================
           DIBUJO
           ==================================================== */

        SDL_SetRenderDrawColor(renderer, 245, 245, 245, 255);
        SDL_RenderClear(renderer);


        /* Trayectoria */

        SDL_SetRenderDrawColor(renderer, 170, 170, 170, 255);

        if (numero_puntos > 1)
            SDL_RenderDrawLines(renderer, trayectoria, numero_puntos);


        /* Brazos */

        SDL_SetRenderDrawColor(renderer, 30, 30, 30, 255);

        for (int i = 0; i < n; i++)
            SDL_RenderDrawLine(
                renderer,
                px[i], py[i],
                px[i + 1], py[i + 1]
            );


        /* Pivote */

        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        dibujar_circulo(renderer, px[0], py[0], 6);


        /* Masas: degradado de azul (primera) a rojo (ultima) */

        for (int i = 1; i <= n; i++)
        {
            double f = (n == 1) ? 1.0 : (double)(i - 1) / (n - 1);

            SDL_SetRenderDrawColor(
                renderer,
                (Uint8)(60  + f * (220 - 60)),
                (Uint8)(110 + f * (70  - 110)),
                (Uint8)(220 + f * (70  - 220)),
                255
            );

            dibujar_circulo(renderer, px[i], py[i], radio_masa);
        }


        SDL_RenderPresent(renderer);


        /* ====================================================
           DATOS EN TERMINAL
           ==================================================== */

        contador++;

        if (contador >= 15)
        {
            contador = 0;

            printf(
                "\nt = %6.2f s | theta1 = %8.2f deg | thetaN = %8.2f deg | E = %12.4f J",
                tiempo,
                estado[0] * 180.0 / PI,
                estado[n - 1] * 180.0 / PI,
                energia_mecanica(estado, &sistema)
            );

            fflush(stdout);
        }


        /* Limitar a ~30 FPS */

        Uint32 ahora = SDL_GetTicks();
        Uint32 transcurrido = ahora - anterior;

        if (transcurrido < 30)
            SDL_Delay(30 - transcurrido);

        anterior = SDL_GetTicks();
    }


    printf("\n\nSimulacion terminada.\n");

    free(trayectoria);
    free(estado);
    free(px);
    free(py);
    sistema_liberar(&sistema);

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
