//
// Created by Hugo Bogaart on 01/10/2024.
//

#ifndef SIMULATION_H
#define SIMULATION_H

#include "Simulation.h"
#include "smath.h"
#include <thread>
#include <mutex>


// func has signature Vector<T, N> -> Vector <T -> N>
template <typename T, size_t N, class Func>
auto RK4 (const Math::Vector<T, N> &in, const T &dt, const Func &func) -> Math::Vector<T, N>
{
        typedef Math::Vector<T, N> Data;
        static_assert(std::is_invocable_r_v<Data, Func, Data>, "Func must have signature Vector<T, N> -> Vector<T, N>");

        const T half_dt = .5 * dt;

        const Data k1 = func(in);
        const Data k2 = func(in + half_dt * k1);
        const Data k3 = func(in + half_dt * k2);
        const Data k4 = func(in + dt * k3);
        return in + dt * (k1 + 2. * (k2 + k3) + k4) / 6.;
}

class DoublePendulumSimulation {

public:
        // phi is the angle with the first mass, psi with te second
        // both relative to the vertical axis
        struct State {
                double phi;
                double psi;
                double phi_deriv;
                double psi_deriv;

        } state;

        // constants
        double l1;
        double l2;

        double m1;
        double m2;


        explicit
        DoublePendulumSimulation (const State &init, double l1, double l2, double m1, double m2)
                : state(init), l1(l1), l2(l2), m1(m1), m2(m2)
        { }


        auto getPhi() const -> double {return state.phi;}
        auto getPsi() const -> double {return state.psi;}
        auto getL1() const -> double {return l1;}
        auto getL2() const -> double {return l2;}


        static constexpr double g = 9.81;


        auto next_it (double dt)
        {
                // y = (phi, psi, u, w)
                // where u := d(phi)/dt and w := d(psi)/dt
                auto diff = [&, this](const Math::Vector<double, 4> &y) -> Math::Vector<double, 4> {
                        const double &m1 = this->m1;
                        const double &m2 = this->m2;
                        const double &l1 = this->l1;
                        const double &l2 = this->l2;

                        const double &phi = y[0];
                        const double &psi = y[1];
                        const double &u   = y[2];
                        const double &w   = y[3];
                        const double phi_psi = phi - psi;
                        const double mu2 = m2 / (m1 + m2);
                        const double sin_phi = std::sin(phi);
                        const double cos_phi = std::cos(phi);
                        const double cos_psi = std::cos(psi);
                        const double sin_phi_psi = std::sin(phi_psi);  // more efficient?
                        const double cos_phi_psi = std::cos(phi_psi);

                        const double u_deriv_up = -m1 * g * sin_phi - m2 * sin_phi_psi * (l1 * u * u * cos_phi_psi + l2 * w * w + g * cos_psi);
                        const double u_deriv = u_deriv_up / (m1 * l1 + m2 * l1 * sin_phi_psi * sin_phi_psi);

                        const double w_deriv_up = l1 * u * u + g * cos_phi + mu2 * l2 * w * w * cos_phi_psi;
                        const double w_deriv   = sin_phi_psi *  w_deriv_up / (l2 - l2 * mu2 * cos_phi_psi * cos_phi_psi);

                        return {u, w, u_deriv, w_deriv};
                };

                const Math::Vector<double, 4> init = {state.phi, state.psi, state.phi_deriv, state.psi_deriv};
                const Math::Vector<double, 4> next = RK4(init, dt, diff);

                state.phi = next[0];
                state.psi = next[1];
                state.phi_deriv = next[2];
                state.psi_deriv = next[3];
        }

        // calculates energy of the system
        // the potential energy is defined to be zero at phi, psi == 0,
        auto getEnergy (const State &state) const -> double
        {
                const double &phi = state.phi;
                const double &psi = state.psi;

                const double &phi_deriv = state.phi_deriv;
                const double &psi_deriv = state.psi_deriv;


                const double cos_phi = std::cos(phi);
                const double cos_psi = std::cos(psi);
                const double cos_phi_psi = std::cos(phi - psi);

                const double V1 = m1 * g * l1 * (1. - cos_phi);
                const double V2 = m2 * g * (l1 + l2 - l1 * cos_phi - l2 * cos_psi);
                const double T1 = 0.5 * m1 * l1 * l1 * phi_deriv * phi_deriv;
                const double T2 = m2 * (0.5 * l1 * l1 * phi_deriv * phi_deriv
                                      + 0.5 * l2 * l2 * psi_deriv * psi_deriv
                                      + l1 * l2 * phi_deriv * psi_deriv * cos_phi_psi);
                return V1 + V2 + T1 + T2;
        }
};

class DPSimulationManager {
public:
        struct DataPoint {
                DoublePendulumSimulation::State state;
                size_t idx;
        };

        double seconds_per_tick;
        DoublePendulumSimulation &sim;
        std::vector<DataPoint> data_points;
        size_t latest_idx = 0;

        DPSimulationManager (DoublePendulumSimulation &simulation, double seconds_per_tick)
                : sim(simulation), seconds_per_tick(seconds_per_tick)
        {
                DataPoint init_dp;
                init_dp.idx = 0;
                init_dp.state = sim.state;
                data_points.push_back(init_dp);
                launch();
        }

        ~DPSimulationManager()
        {
                {
                        std::lock_guard _(mtx);
                        kill_thread = true;
                }
                sim_runner.join();
        }

        // sets the number of datapoints we aim to have in the vector
        void set_target_size (size_t n)
        {
                std::lock_guard _(this->mtx);
                target_size = n;
        }

        DataPoint extract_data_point ()
        {
                // first we need to make sure there is a datapoint
                while (true) {
                        bool has_dp;
                        {
                                std::lock_guard _(this->mtx);
                                has_dp = !data_points.empty();
                        }
                        if (has_dp)
                                break;

                        // we wait for thread to finish the tick
                        // std::this_thread::sleep_for(std::chrono::milliseconds(10));
                        std::this_thread::yield();
                }
                std::lock_guard _(this->mtx);
                DataPoint dp = data_points.front();
                data_points.erase(data_points.begin());
                return dp;
        }

private:

        void launch ()
        {
                auto thread_loop = [this]() {

                        while (true) {
                                bool kill;
                                bool target_satisfied;
                                {
                                        std::lock_guard _(this->mtx);
                                        kill = this->kill_thread;
                                        target_satisfied = data_points.size() >= target_size;
                                }
                                if (kill)
                                        return;

                                if (target_satisfied) {
                                        // std::this_thread::sleep_for(std::chrono::seconds(1));
                                        std::this_thread::yield();
                                } else {
                                        // we run the sim and add it to the target
                                        sim.next_it(seconds_per_tick);
                                        ++latest_idx;
                                        DataPoint dp;
                                        dp.idx = latest_idx;
                                        dp.state = sim.state;
                                        {
                                                std::lock_guard _(this->mtx);
                                                data_points.push_back(dp);
                                        }
                                }
                        }
                };
                sim_runner = std::thread(thread_loop);
        }


        std::mutex mtx;
        std::thread sim_runner;
        size_t target_size = 10;
        bool kill_thread = false;
};



#endif //SIMULATION_H
