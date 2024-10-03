#include <iostream>

extern "C" {
#include "raylib.h"
#include "rlgl.h"
#include "raymath.h"
}

#include <chrono>
#include <random>
#include <optional>
#include <vector>
#include <array>

#include "smath.h"
#include "Simulation.h"

using Math::Vector;

auto operator+ (const Vector2& lhs, const Vector2& rhs) -> Vector2
{
        return Vector2{lhs.x + rhs.x, lhs.y + rhs.y};
}

auto operator- (const Vector2& lhs, const Vector2& rhs) -> Vector2
{
        return Vector2{lhs.x - rhs.x, lhs.y - rhs.y};
}
auto operator* (const float& l, const Vector2& rhs) -> Vector2
{
        return Vector2{l * rhs.x, l * rhs.y};
}
auto operator* (const Vector2& lhs, const float& l) -> Vector2
{
        return l * lhs;
}

auto operator/ (const Vector2& lhs, const float& l) -> Vector2
{
        return Vector2{.x = lhs.x / l, .y = lhs.y / l};
}

auto drawStaticGrid (const Camera2D &camera, float spacing_px)
{
        const Color color = GRAY;
        const float spacing = spacing_px / camera.zoom;
        const int pix_x = GetScreenWidth();
        const int pix_y = GetScreenHeight();
        Vector2 wdw = GetScreenToWorld2D({0, 0}, camera);
        const float min_x = wdw.x;
        const float max_x = min_x + spacing * (float)pix_x;
        const float min_y = wdw.y;
        const float max_y = min_y + spacing * (float)pix_y;
        const float center_x = (min_x + max_x) / 2;
        const float center_y = (min_y + max_y) / 2;
        for (float x = center_x; x > min_x; x -= spacing) {
                DrawLineV({x, min_y}, {x, max_y}, color);
        }
        for (float x = center_x; x < max_x; x += spacing) {
                DrawLineV({x, min_y}, {x, max_y}, color);
        }
        for (float y = center_y; y > min_y; y -= spacing) {
                DrawLineV({min_x, y}, {max_x, y}, color);
        }
        for (float y = center_y; y < max_y; y += spacing) {
                DrawLineV({min_x, y}, {max_x, y}, color);
        }
}

struct GuiState {
        bool show_grid = true;
        bool show_energy = false;
        Color background_color = RAYWHITE;
        bool trace = false;
        bool show_lines = true;
        Color line_color = BLACK;
};

auto flip (bool &b) {b = !b;}


int main()
{
        constexpr int target_fps = 60;
        constexpr int target_tps = 480;
        constexpr double target_spt = 1. / static_cast<double>(target_tps);
        constexpr int tpf = target_tps / target_fps;


        // PendulumSimulation sim(.7, 0, 3, 1, false);
        DoublePendulumSimulation::State init_state{};
        init_state.phi = 0; //.6;
        init_state.psi = M_PI; //2;
        init_state.phi_deriv = 2.1; //0;
        init_state.psi_deriv = 0;
        float m1_metric = 10;
        float m2_metric = 6.5;
        float l1_metric = 2;
        float l2_metric = 1.3;

        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution dis(-.2, .2);
        double random = dis(gen);
        init_state.psi_deriv += random;

        DoublePendulumSimulation dsim(init_state, l1_metric, l2_metric, m1_metric, m2_metric);
        DPSimulationManager sim_manager(dsim, target_spt);
        sim_manager.set_target_size(target_tps * 10);

        // return 0;

        // Initialization
        //--------------------------------------------------------------------------------------
        const int screenWidth = 800;
        const int screenHeight = 450;
        const double units_per_meter = 100;


        SetConfigFlags(FLAG_WINDOW_RESIZABLE); // Window configuration flags
        InitWindow(screenWidth, screenHeight, "Simulatie");

        Camera2D camera = {};
        camera.zoom = 1.0f;

        int zoomMode = 0; // 0-Mouse Wheel, 1-Mouse Move


        SetTargetFPS(target_fps); // Set our game to run at 60 frames-per-second
        //--------------------------------------------------------------------------------------

        GuiState gstate;

        struct TraceData {
                Vector2 pos;
                int value;
                float rad;
                int ball;
        };
        std::vector<TraceData> trace1, trace2, trace1w, trace2w;
        trace1.reserve(1000 * 100);
        trace2.reserve(1000 * 100);
        trace1w.reserve(1000 * 100);
        trace2w.reserve(1000 * 100);

        const float rect_height_units = 50;
        const float rect_width_units = 120;

        // Main game loop
        while (!WindowShouldClose()) // Detect window close button or ESC key
        {
                // CAMERA ::
                // offset is the pixel in the screen that is considered the "screen origin"
                // target is the world coordinate that is mapped to the screen origin

                // Update
                //----------------------------------------------------------------------------------

                camera.offset = {(float)GetScreenWidth() / 2, rect_height_units};
                camera.target = {0, 0};
                Vector2 rect_attachment_point = {0, 0};

                if (zoomMode == 0) {
                        // Zoom based on mouse wheel
                        float wheel = GetMouseWheelMove();
                        if (wheel != 0) {
                                // Get the world point that is under the mouse
                                // Vector2 mouseWorldPos = GetScreenToWorld2D(GetMousePosition(), camera);

                                // Set the offset to where the mouse is
                                // camera.offset = GetMousePosition();

                                // Set the target to match, so that the camera maps the world space point
                                // under the cursor to the screen space point under the cursor at any zoom
                                // camera.target = mouseWorldPos;

                                // Zoom increment
                                float scaleFactor = 1.0f + (0.25f * fabsf(wheel));
                                if (wheel < 0)
                                        scaleFactor = 1.0f / scaleFactor;
                                camera.zoom = Clamp(camera.zoom * scaleFactor, 0.125f, 64.0f);
                        }
                }

                // key events
                if (IsKeyPressed(KEY_E)) {
                        flip(gstate.show_energy);
                }

                if (IsKeyPressed(KEY_G)) {
                        flip(gstate.show_grid);
                }
                if (IsKeyPressed(KEY_T)) {
                        flip(gstate.trace);
                        if (gstate.trace) {
                                gstate.background_color = BLACK;
                                gstate.line_color = {255, 255, 255, 100};
                        } else {
                                gstate.background_color = RAYWHITE;
                                gstate.line_color = BLACK;
                                trace1.clear();
                                trace2.clear();
                                trace1w.clear();
                                trace2w.clear();
                        }
                }
                if (IsKeyPressed(KEY_L)) {
                        flip(gstate.show_lines);
                }

                // some coordinates we use


                // black rectangle

                Vector2 rectV = rect_attachment_point + Vector2{.x = -rect_width_units / 2, .y = -rect_height_units};

                std::array<DPSimulationManager::DataPoint, tpf> sim_data{};

                for (int i = 0; i < tpf; i++) {
                        sim_data[i] = sim_manager.extract_data_point();
                }
                const auto &current_state = sim_data.back().state;

                const float l1 = dsim.getL1() * units_per_meter;
                const float l2 = dsim.getL2() * units_per_meter;

                // world
                std::array<Vector2, tpf> positions1{};
                std::array<Vector2, tpf> positions2{};

                for (int i = 0; i < tpf; i++) {
                        const float phi = sim_data[i].state.phi;
                        const float psi = sim_data[i].state.psi;
                        positions1[i] = rect_attachment_point + Vector2{l1 * std::sin(phi), l1 * std::cos(phi)};
                        positions2[i] = positions1[i] + Vector2{l2 * std::sin(psi), l2 * std::cos(psi)};
                }

                Vector2 ball1V = positions1.back();
                Vector2 ball2V = positions2.back();

                // radii of the balls in world units
                const float rad1  = 4;
                const float rad2  = 6;

                const float rad1_trace  = 2;
                const float rad2_trace  = 4;

                // if we keep track of the trace, we decrease all values by one
                // and make the new point
                constexpr int trace_lifetime1 = 40;
                constexpr int trace_lifetime2 = 300;

                constexpr float rad_decay1 = 0.9; // 0.9865;
                constexpr float rad_decay2 = 0.9865; // 0.9865;


                if (gstate.trace) {
                        // all the traces we keep go into trace<i>w
                        trace1w.clear();
                        trace2w.clear();
                        for (auto &s : trace1) {
                                --s.value;
                                s.rad *= rad_decay1;
                                if (s.value > 0) {
                                        trace1w.push_back(s);
                                }
                        }
                        for (auto &s : trace2) {
                                --s.value;
                                s.rad *= rad_decay2;
                                if (s.value > 0) {
                                        trace2w.push_back(s);
                                }
                        }
                        std::swap(trace1, trace1w);
                        std::swap(trace2, trace2w);

                        // so now all the good ones are in trace<i>w

                        // now we add all the reached positions to the traces
                        // we can also interpolate between the datapoints for more balls

                        constexpr size_t num_interpolate = 3;
                        auto interpolate = [&](const TraceData &td1, const TraceData &td2) -> std::array<TraceData, num_interpolate> {
                                std::array<TraceData, num_interpolate> result;
                                Vector2 segment = (td2.pos - td1.pos) / static_cast<float>(num_interpolate + 1);
                                float dr = (td2.rad - td1.rad) / static_cast<float>(num_interpolate + 1);
                                for (size_t i = 0; i < num_interpolate; i++) {
                                        result[i] = td2;
                                        result[i].pos = td1.pos + static_cast<float>(i + 1) * segment;
                                        result[i].rad = td1.rad + static_cast<float>(i + 1) * dr;
                                }
                                return result;
                        };


                        for (int i = 0; i < tpf; i++) {
                                TraceData td1 = {positions1[i], trace_lifetime1, rad1_trace, 0};
                                TraceData td2 = {positions2[i], trace_lifetime2, rad2_trace, 1};

                                std::optional<std::array<TraceData, num_interpolate>> data1{}, data2{};
                                if (!trace1.empty())
                                        data1 = interpolate(trace1.back(), td1);
                                if (!trace2.empty())
                                        data2 = interpolate(trace2.back(), td2);

                                if (data1) {
                                        for (const auto &tr : *data1) {
                                                trace1.push_back(tr);
                                        }
                                }

                                if (data2) {
                                        for (const auto &tr : *data2) {
                                                trace2.push_back(tr);
                                        }
                                }

                                trace1.push_back(td1);
                                trace2.push_back(td2);
                        }
                }


                //----------------------------------------------------------------------------------

                // Draw
                //----------------------------------------------------------------------------------
                BeginDrawing();
                ClearBackground(gstate.background_color);
                BeginMode2D(camera);

                if (gstate.show_lines) {
                        DrawLineEx(rect_attachment_point, ball1V, 1., gstate.line_color);
                        DrawLineEx(ball1V, ball2V, 1., gstate.line_color);
                }

                // grid is ignored if we trace
                if (gstate.trace) {
                        DrawCircleV(rect_attachment_point, .5, WHITE);

                        for (auto &s : trace1) {
                                float frac = (float) s.value / (float) trace_lifetime1;
                                unsigned char rgb = 255 * frac;
                                Color col = {170, 210, 255, rgb};
                                DrawCircleV(s.pos, s.rad, col);
                        }
                        for (auto &s : trace2) {
                                float frac = (float) s.value / (float) trace_lifetime2;
                                unsigned char rgb = 255 * frac * .4;
                                Color col = Color{255, 255, 255, rgb};
                                DrawCircleV(s.pos, s.rad, col);
                        }
                        
                } else {
                        if (gstate.show_grid)
                                drawStaticGrid(camera, 30);


                        DrawRectangleV(rectV, {.x = rect_width_units, .y = rect_height_units}, BLACK);

                        const Color ball_col = GRAY;

                        DrawCircleV(rect_attachment_point, 2, BLACK);
                        DrawCircleV(ball1V, rad1, ball_col);
                        DrawCircleV(ball2V, rad2, ball_col);
                }
                EndMode2D();

                if (gstate.show_energy) {
                        double energy = dsim.getEnergy(current_state);
                        std::string txt = "E: ";
                        txt.append(std::to_string(energy));
                        txt.append(" J");
                        DrawText(txt.c_str(), 10, 10, 22, GREEN);

                        std::string fps_str = "FPS: ";
                        fps_str += std::to_string(GetFPS());
                        DrawText(fps_str.c_str(), 10, 40, 22, GREEN);

                        std::string tps_str = "TPS: ";
                        tps_str += std::to_string(target_tps);
                        DrawText(tps_str.c_str(), 10, 70, 22, GREEN);

                        std::string sz_str = "CACHE: ";
                        sz_str += std::to_string(sim_manager.data_points.size());
                        DrawText(sz_str.c_str(), 10, 100, 22, GREEN);
                }


                EndDrawing();
                //----------------------------------------------------------------------------------
        }

        // De-Initialization
        //--------------------------------------------------------------------------------------
        CloseWindow(); // Close window and OpenGL context
        //--------------------------------------------------------------------------------------
        return 0;
}