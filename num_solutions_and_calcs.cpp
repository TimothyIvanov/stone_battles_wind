#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <tuple>
#include <iomanip>
#include <sstream>

// Function to calculate speed relative to wind
double vr(double vx, double vy, double vw)
{
    return sqrt(pow(vx - vw, 2) + pow(vy, 2));
}

// Function to calculate x-component of the net force
double Fx(double vx, double vy, double vw, double alpha, double m)
{
    return -alpha * vr(vx, vy, vw) * (vx - vw) / m;
}

// Function to calculate y-component of the net force
double Fy(double vx, double vy, double vw, double alpha, double m, double g)
{
    return -g - alpha * vr(vx, vy, vw) * vy / m;
}

// Function for a single step in solving ODEs
// (Adding a structure for ODEstep interferes with solveODE for some reason, not enough time to figure out why)
std::tuple<double, double, double> ODEstep(double t, double vx, double vy, double vw, double alpha, double m, double g, double dt)
{
    double vx0 = vx, vy0 = vy;
    t = t + dt;
    vx = vx0 + dt * Fx(vx0, vy0, vw, alpha, m);
    vy = vy0 + dt * Fy(vx0, vy0, vw, alpha, m, g);
    return std::make_tuple(t, vx, vy);
}

// Structure to store ODE solution
struct ODEResult {
    std::vector<double> t_list;
    std::vector<double> vx_list;
    std::vector<double> vy_list;
};

// Function to solve ODE and store the results
ODEResult solveODE(double v0, double th0, double vw, double alpha, double m, double g, double tmax, double dt)
{
    double t = 0.0;
    double vx = v0 * std::cos(M_PI / 180.0 * th0);
    double vy = v0 * std::sin(M_PI / 180.0 * th0);

    ODEResult result;

    // Iterate through time steps
    while (t <= tmax)
    {
        result.t_list.push_back(t);
        result.vx_list.push_back(vx);
        result.vy_list.push_back(vy);

        std::tie(t, vx, vy) = ODEstep(t, vx, vy, vw, alpha, m, g, dt);
    }
    return result;
}

// Structure to store position coordinates
struct PositionResult {
    std::vector<double> x;
    std::vector<double> y;
};

// Function to calculate positions from velocities
PositionResult calcPos(const double &x0, const double &y0, const std::vector<double> &t, const std::vector<double> &vx, const std::vector<double> &vy)
{
    std::vector<double> x(t.size());
    std::vector<double> y(t.size());

    x[0] = x0;
    y[0] = y0;

    // Integrate velocities to get positions
    for (size_t i = 1; i < t.size(); ++i)
    {
        double dt = t[i] - t[i - 1];
        x[i] = x[i - 1] + vx[i] * dt;
        y[i] = y[i - 1] + vy[i] * dt;
    }

    PositionResult result;
    result.x = x;
    result.y = y;

    return result;
}

// Function to save simulation data to files
void saveData(const std::string &filename, const std::string &title, const std::vector<double> &t, const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &vx, const std::vector<double> &vy)
{
    std::ofstream file(filename);
    if (file.is_open())
    {
        file << "#" << title << "\n";
        for (size_t i = 0; i < t.size(); ++i)
        {
            file << std::fixed << std::setprecision(6) << t[i] << " " << x[i] << " " << y[i] << " " << vx[i] << " " << vy[i] << "\n";
        }
        file.close();
    }
    else
    {
        std::cerr << "Error opening file: " << filename << std::endl;
    }
}

// Function to format & convert double values because C++ likes to make me suffer
std::string formatDouble(double value, int precision)
{
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(precision) << value;
    
    return ss.str();
}

// Main simulation function
// I'm aware it's a little repetitive and I can write a function to make it more concise
// No time, however, sorry :)
int main()
{
    // Wind velocities and time steps/dt's
    std::vector<double> wind_velocities = {2.1, 0, -2.1};
    std::vector<double> time_steps = {0.01, 0.005};

    // Environment variables
    const double alpha = 2.5; // air drag strength in kg/m
    const double m = 2.8;     // mass in kg
    const double g = 9.81;    // acceleration of gravity in m/s^2

    // Initial conditions
    const double x0 = 0.0;    // initial horizontal coordinate in m
    const double h0 = 2.1;    // initial vertical coordinate in m
    const double v0 = 6.9;    // initial speed in m/s
    const double th0 = 34;    // initial direction in degrees

    const double maxt = 5;    // maximum time delta

    // Loop over wind velocities and time steps
    for (auto vw : wind_velocities)
    {
        for (auto dt : time_steps)
        {
            // Get DE solutions
            ODEResult odeResult = solveODE(v0, th0, vw, alpha, m, g, maxt, dt);

            // Format wind and time step
            std::string vw_str = formatDouble(vw, 1);
            std::string dt_str = formatDouble(dt, 3);

            // Get object trajectory
            PositionResult posResult = calcPos(x0, h0, odeResult.t_list, odeResult.vx_list, odeResult.vy_list);

            // Save formatted data to files
            saveData("vw_" + vw_str + "_dt_" + dt_str + ".dat", "t       x        y        vx        vy", odeResult.t_list, posResult.x, posResult.y, odeResult.vx_list, odeResult.vy_list);
        }
    }

    // Initial conditions and parameters for (0,0) landing trial
    double t_h0 = 7.7;
    double t_dt = 0.005;
    std::vector<double> t_wind_velocities = {-1, -0.5, -0.7, -0.68, -0.6775};

    //Loop over trial wind velocities
    for (auto t_vw : t_wind_velocities)
    {
        // Get DE solutions
        ODEResult trial_odeResult = solveODE(v0, th0, t_vw, alpha, m, g, maxt, t_dt);

        // Format wind and time step
        std::string t_vw_str = formatDouble(t_vw, 4);
        std::string t_dt_str = formatDouble(t_dt, 3);

        // Get object trajectory
        PositionResult posResult = calcPos(x0, t_h0, trial_odeResult.t_list, trial_odeResult.vx_list, trial_odeResult.vy_list);

        // Save formatted data to files
        saveData("t_vw_" + t_vw_str + "_t_dt_" + t_dt_str + ".dat", "t       x        y        vx        vy", trial_odeResult.t_list, posResult.x, posResult.y, trial_odeResult.vx_list, trial_odeResult.vy_list);
    }
}