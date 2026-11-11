#include "raylib.h"
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define MAX_TRAIL 1200
#define SUBSTEPS 12

// --- Estructura del Péndulo Doble 3D ---
typedef struct {
    Vector3 r1, r1_old; // Posición de la masa 1 (actual y previa)
    Vector3 r2, r2_old; // Posición de la masa 2 (actual y previa)
    float m1, m2;       // Masas
    float L1, L2;       // Longitudes de las varillas
    float damping;      // Resistencia del aire (permite ver la formación de ondas)
} PenduloDoble3D;

// --- Buffer circular para el rastro (trail) ---
typedef struct {
    Vector3 puntos[MAX_TRAIL];
    int contador;
    int indice;
} Rastro;

// --- Operaciones Vectoriales Básicas ---
static inline Vector3 VAdd(Vector3 a, Vector3 b) { return (Vector3){ a.x + b.x, a.y + b.y, a.z + b.z }; }
static inline Vector3 VSub(Vector3 a, Vector3 b) { return (Vector3){ a.x - b.x, a.y - b.y, a.z - b.z }; }
static inline Vector3 VScale(Vector3 a, float s)  { return (Vector3){ a.x * s, a.y * s, a.z * s }; }
static inline float VDot(Vector3 a, Vector3 b)    { return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline float VLengthSq(Vector3 a)          { return VDot(a, a); }
static inline Vector3 VLerp(Vector3 a, Vector3 b, float t) {
    return (Vector3){ a.x + (b.x - a.x)*t, a.y + (b.y - a.y)*t, a.z + (b.z - a.z)*t };
}

// Inicializa las condiciones del péndulo
void ReiniciarPendulo(PenduloDoble3D *p, Rastro *trail) {
    p->m1 = 1.0f;
    p->m2 = 1.0f;
    p->L1 = 2.4f;
    p->L2 = 2.0f;
    p->damping = 0.012f; // Fricción suave para apreciar ondas concéntricas

    // Ángulos iniciales desviados en 3D (para generar movimiento en X, Y y Z)
    float theta1 = 1.45f, phi1 = 0.35f;
    p->r1 = (Vector3){
        p->L1 * sinf(theta1) * cosf(phi1),
        -p->L1 * cosf(theta1),
        p->L1 * sinf(theta1) * sinf(phi1)
    };

    float theta2 = 1.20f, phi2 = -0.70f;
    Vector3 rel2 = {
        p->L2 * sinf(theta2) * cosf(phi2),
        -p->L2 * cosf(theta2),
        p->L2 * sinf(theta2) * sinf(phi2)
    };
    p->r2 = VAdd(p->r1, rel2);

    // Velocidades iniciales con impulso lateral
    Vector3 v1 = { 0.0f, 0.0f, 1.2f };
    Vector3 v2 = { 0.4f, 0.0f, -1.0f };

    float dt = (1.0f / 60.0f) / SUBSTEPS;
    Vector3 g = { 0.0f, -9.81f, 0.0f };
    p->r1_old = VAdd(VSub(p->r1, VScale(v1, dt)), VScale(g, 0.5f * dt * dt));
    p->r2_old = VAdd(VSub(p->r2, VScale(v2, dt)), VScale(g, 0.5f * dt * dt));

    trail->contador = 0;
    trail->indice = 0;
}

// Integración Verlet con algoritmo SHAKE (conserva rigidez y energía)
void ActualizarFisica(PenduloDoble3D *p, float dt) {
    Vector3 g = { 0.0f, -9.81f, 0.0f };

    // Velocidades estimadas para la amortiguación por aire
    Vector3 v1 = VScale(VSub(p->r1, p->r1_old), 1.0f / dt);
    Vector3 v2 = VScale(VSub(p->r2, p->r2_old), 1.0f / dt);

    Vector3 a1 = VSub(g, VScale(v1, p->damping));
    Vector3 a2 = VSub(g, VScale(v2, p->damping));

    // Paso de Verlet no restringido
    Vector3 s1 = VAdd(VSub(VScale(p->r1, 2.0f), p->r1_old), VScale(a1, dt * dt));
    Vector3 s2 = VAdd(VSub(VScale(p->r2, 2.0f), p->r2_old), VScale(a2, dt * dt));

    Vector3 r1_new = s1;
    Vector3 r2_new = s2;

    // Iteraciones SHAKE para satisfacer la distancia de las varillas
    for (int it = 0; it < 20; it++) {
        // Restricción 1: |r1| == L1 (Pivote fijo en el origen 0,0,0)
        float diff1 = VLengthSq(r1_new) - p->L1 * p->L1;
        if (fabsf(diff1) > 1e-7f) {
            float delta1 = -diff1 / (2.0f * VDot(r1_new, p->r1));
            r1_new = VAdd(r1_new, VScale(p->r1, delta1 / p->m1));
        }

        // Restricción 2: |r2 - r1| == L2
        Vector3 d_new = VSub(r2_new, r1_new);
        Vector3 d_old = VSub(p->r2, p->r1);
        float diff2 = VLengthSq(d_new) - p->L2 * p->L2;
        if (fabsf(diff2) > 1e-7f) {
            float denom = 2.0f * VDot(d_new, d_old) * (1.0f / p->m1 + 1.0f / p->m2);
            float delta2 = -diff2 / denom;
            r1_new = VSub(r1_new, VScale(d_old, delta2 / p->m1));
            r2_new = VAdd(r2_new, VScale(d_old, delta2 / p->m2));
        }
    }

    p->r1_old = p->r1;
    p->r1 = r1_new;
    p->r2_old = p->r2;
    p->r2 = r2_new;
}

// Agrega la posición actual al rastro
void AgregarRastro(Rastro *t, Vector3 pos) {
    t->puntos[t->indice] = pos;
    t->indice = (t->indice + 1) % MAX_TRAIL;
    if (t->contador < MAX_TRAIL) t->contador++;
}

int main(void) {
    const int anchoVentana = 1280;
    const int altoVentana = 720;
    InitWindow(anchoVentana, altoVentana, "Pendulo Doble 3D: Ondas (X-Z) y Perspectiva de Profundidad (Y-Z)");

    PenduloDoble3D pendulo;
    Rastro rastro;
    ReiniciarPendulo(&pendulo, &rastro);

    // --- Definición de Cámaras ---
    // Perspectiva 0: Frontal (profundidad en Z y oscilación en Y)
    Vector3 camPosFrontal   = { 0.0f, -2.0f, 7.5f };
    Vector3 camTgtFrontal   = { 0.0f, -2.2f, 0.0f };
    Vector3 camUpFrontal    = { 0.0f, 1.0f, 0.0f };

    // Perspectiva 1: Cenital / Superior (mirando las ondas en el plano X-Z)
    Vector3 camPosCenital   = { 0.0f, 8.8f, 0.001f }; // Pequeño offset en Z para evitar gimbal con 'up'
    Vector3 camTgtCenital   = { 0.0f, 0.0f, 0.0f };
    Vector3 camUpCenital    = { 0.0f, 0.0f, -1.0f };

    int modoVista = 0; // 0 = Frontal (Profundidad Z), 1 = Cenital (Ondas X-Z)

    Camera3D camara = { 0 };
    camara.position = camPosFrontal;
    camara.target   = camTgtFrontal;
    camara.up       = camUpFrontal;
    camara.fovy     = 55.0f;
    camara.projection = CAMERA_PERSPECTIVE;

    SetTargetFPS(60);

    while (!WindowShouldClose()) {
        // --- Controles de Usuario ---
        if (IsKeyPressed(KEY_SPACE) || IsKeyPressed(KEY_TAB)) {
            modoVista = 1 - modoVista; // Alterna la perspectiva en tiempo real
        }
        if (IsKeyPressed(KEY_R)) {
            ReiniciarPendulo(&pendulo, &rastro);
        }
        if (IsKeyPressed(KEY_D)) {
            pendulo.damping = (pendulo.damping > 0.0f) ? 0.0f : 0.012f;
        }

        // --- Transición Suave de Cámara (Interpolación) ---
        Vector3 targetPos = (modoVista == 0) ? camPosFrontal : camPosCenital;
        Vector3 targetTgt = (modoVista == 0) ? camTgtFrontal : camTgtCenital;
        Vector3 targetUp  = (modoVista == 0) ? camUpFrontal  : camUpCenital;

        camara.position = VLerp(camara.position, targetPos, 0.08f);
        camara.target   = VLerp(camara.target, targetTgt, 0.08f);
        camara.up       = VLerp(camara.up, targetUp, 0.08f);

        // --- Actualización de Física por Subpasos ---
        float dtFrame = GetFrameTime();
        if (dtFrame > 0.033f) dtFrame = 0.033f; // Previene saltos bruscos
        float subDt = dtFrame / SUBSTEPS;

        for (int i = 0; i < SUBSTEPS; i++) {
            ActualizarFisica(&pendulo, subDt);
        }
        AgregarRastro(&rastro, pendulo.r2);

        // --- Renderizado Gráfico ---
        BeginDrawing();
        ClearBackground((Color){ 16, 18, 26, 255 }); // Azul oscuro profundo

        BeginMode3D(camara);

            // Rejilla de referencia en el piso
            DrawGrid(16, 1.0f);

            // 1. Dibujar el Rastro (Ondas tridimensionales y proyección en el suelo)
            float nivelPiso = -4.8f;
            for (int i = 0; i < rastro.contador - 1; i++) {
                int idx1 = (rastro.indice - rastro.contador + i + MAX_TRAIL) % MAX_TRAIL;
                int idx2 = (idx1 + 1) % MAX_TRAIL;

                float factorEdad = (float)i / (float)rastro.contador;
                Color col3D = ColorAlpha(SKYBLUE, factorEdad * 0.85f);
                Color colPiso = ColorAlpha(MAGENTA, factorEdad * 0.45f);

                // Curva 3D en el aire
                DrawLine3D(rastro.puntos[idx1], rastro.puntos[idx2], col3D);

                // Proyección de la "onda de arena" en el piso (ideal para vista cenital)
                Vector3 pFloor1 = { rastro.puntos[idx1].x, nivelPiso, rastro.puntos[idx1].z };
                Vector3 pFloor2 = { rastro.puntos[idx2].x, nivelPiso, rastro.puntos[idx2].z };
                DrawLine3D(pFloor1, pFloor2, colPiso);
            }

            // 2. Dibujar el Péndulo
            Vector3 pivote = { 0.0f, 0.0f, 0.0f };

            // Varillas (Cilindros 3D con grosor uniforme en el espacio)
            DrawCylinderEx(pivote, pendulo.r1, 0.035f, 0.035f, 8, LIGHTGRAY);
            DrawCylinderEx(pendulo.r1, pendulo.r2, 0.035f, 0.035f, 8, GRAY);

            // Esferas (Masa 1 y Masa 2)
            DrawSphere(pivote, 0.09f, DARKGRAY);       // Soporte pivote
            DrawSphere(pendulo.r1, 0.18f, ORANGE);     // Masa 1
            DrawSphere(pendulo.r2, 0.24f, GOLD);       // Masa 2 (Extremo)

            // Línea vertical tenue para guiar la profundidad hacia el piso
            DrawLine3D(pendulo.r2, (Vector3){ pendulo.r2.x, nivelPiso, pendulo.r2.z }, ColorAlpha(WHITE, 0.2f));

        EndMode3D();

        // --- Interfaz de Usuario 2D (HUD) ---
        DrawRectangle(15, 15, 410, 195, ColorAlpha(BLACK, 0.75f));
        DrawRectangleLines(15, 15, 410, 195, DARKGRAY);

        DrawText("PÉNDULO DOBLE 3D: DINÁMICA Y PERSPECTIVA", 25, 25, 16, RAYWHITE);

        if (modoVista == 0) {
            DrawText("Modo Actual: VISTA FRONTAL (Profundidad Z)", 25, 50, 15, YELLOW);
            DrawText("• Observa el eje Y (arriba/abajo)", 25, 70, 13, LIGHTGRAY);
            DrawText("• Al viajar en Z, el péndulo se AGRANDA al acercarse", 25, 88, 13, GREEN);
            DrawText("  y se ENCOGE al alejarse en el fondo.", 25, 104, 13, GREEN);
        } else {
            DrawText("Modo Actual: VISTA CENITAL (Ondas X-Z)", 25, 50, 15, MAGENTA);
            DrawText("• Mirando desde arriba hacia el plano horizontal", 25, 70, 13, LIGHTGRAY);
            DrawText("• Se aprecian las curvas y ondas compuestas", 25, 88, 13, SKYBLUE);
            DrawText("  similares a figuras de Lissajous / péndulo de arena.", 25, 104, 13, SKYBLUE);
        }

        // Datos en tiempo real
        DrawText(TextFormat("Posición Extremo -> Y: %.2f | Z (Profundidad): %.2f", pendulo.r2.y, pendulo.r2.z), 25, 130, 14, RAYWHITE);
        DrawText(TextFormat("Fricción Aire: %s", (pendulo.damping > 0.0f ? "ACTIVADA (Forma Ondas)" : "DESACTIVADA (Caos Puro)"), 25, 150, 13, ORANGE);
        DrawText("[ESPACIO]: Cambiar Vista | [R]: Reiniciar | [D]: Toggle Fricción", 25, 175, 13, WHITE);

        EndDrawing();
    }

    CloseWindow();
    return 0;
}
