#include <random>
#include <cmath>
#include <iostream>

void monte_carlo();
void print_substrate();
void analyze_surface();
int count_op(int x, int y);
int count_ip(int x, int y);
void move_op(int x, int y);
void select(int x, int y);
void select_move(int x, int y, int sel);
void add_nn(int flux);
void add_elem(int x, int y);
void rem_elem(int x, int y);

struct nn
{
    uint32_t height = 0;
    bool occupied = false;
    bool vacant = true;
};

main()
{
    using namespace std;
    auto x = 0;
    //std::cout << "hello world: " << std::endl;

    monte_carlo();
    print_substrate();
    analyze_surface();
    return x;
}

//constants
static constexpr float k_b = 8.617333e-5f;    //eV-s
static constexpr float h_bar = 6.582911e-15f; // ev/K

//sim params
static constexpr int iter_max = 100; // iterations
static constexpr int x_max = 25;     //units
static constexpr int y_max = 25;     //units
static constexpr int flux = 100;

//substrate object
nn plant[x_max][y_max];

void monte_carlo()
{
    //basic calculations

    float temp = 900.0f; //Kelvin
    float tau = .0001;

    //base rates
    float dep_rate = 1.0f;                       // in ML/s
    float base_rate = 2.0f * k_b * temp / h_bar; // in Hz

    //base energies
    float E_diff_ip = 1.75f;  // in-plane energies
    float E_diff_op = 0.4f;   // in-plane energies
    float E_desorb_ip = 2.0f; // out-of-plane energies
    float E_desorb_op = 0.2f; // out-of-plane energies

    //equations
    float diff_rate = base_rate * exp(-E_diff_ip / (k_b * temp));
    float desorb_rate = base_rate * exp(-E_desorb_ip / (k_b * temp));

    float diff_rate_max = base_rate * exp(-E_diff_ip * 4.0f / (k_b * temp));
    float desorb_rate_max = base_rate * exp(-E_desorb_ip * 4.0f / (k_b * temp));

    float diff_prob = diff_rate * tau;
    float desorb_prob = desorb_rate * tau;

    float diff_prob_max = diff_rate_max * tau;
    float desorb_prob_max = desorb_rate_max * tau;

    std::cout << "baseline diff rate: " << diff_rate << std::endl;
    std::cout << "baseline desorb rate: " << desorb_rate << std::endl;

    std::cout << "baseline diff P: " << diff_prob << std::endl;
    std::cout << "baseline desorb P: " << desorb_prob << std::endl;

    std::cout << "max diff P: " << diff_prob_max << std::endl;
    std::cout << "max desorb P: " << desorb_prob_max << std::endl;

    //    nn plant[x_max][y_max];
    for (int idx = 0; idx < iter_max; idx++)
    {
        add_nn(flux);
        for (int idx_x = 1; idx_x < x_max - 1; idx_x++)
        {
            for (int idx_y = 1; idx_y < y_max - 1; idx_y++)
            {
                if (plant[idx_x][idx_y].occupied == true)
                {
                    float r = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
                    diff_prob = base_rate * exp(-E_diff_ip * static_cast<float>(count_ip(idx_x, idx_y)) / (k_b * temp));
                    if (r > diff_prob)
                    {
                        std::cout << "start iter: " << idx << ", x: " << idx_x << ", y: " << idx_y << ", P: " << r << std::endl;
                        select(idx_x, idx_y); //atom leaves site and moves to a new one
                        std::cout << "end iter: " << idx << ", x: " << idx_x << ", y: " << idx_y << ", P: " << r << std::endl;
                        //find which new site to locate to
                    }
                    if (r > desorb_prob)
                    {
                        //select(idx_x, idx_y); //adatom leaves lattice
                    }
                }
            }
        }
    }
}

//in-plane nearest-neighbor selection
int count_ip(int x, int y)
{
    //std::cout << "(" << x << "," << y << ") ";
    int count = 0;

    if (plant[x][y - 1].occupied == true)
    {
        //    std::cout << "A, ";
        count++;
    }
    if (plant[x + 1][y].occupied == true)
    {
        //    std::cout << "B, ";
        count++;
    }
    if (plant[x][y + 1].occupied == true)
    {
        //   std::cout << "C, ";
        count++;
    }
    if (plant[x - 1][y].occupied == true)
    {
        //   std::cout << "D ";
        count++;
    }
    //std::cout << "are occupied" << std::endl;

    return count;
}

//out-of-plane nearest-neighbor selection
int count_op(int x, int y)
{
    int count = 0;
    if (plant[x][y].vacant == true)
    {
        count++;
    }
    return count;
}

void select(int x, int y)
{
    int max_ip = 4;
    int max_op = 1;

    int num_options_tot = max_ip + count_op(x, y);

    int select = rand() % num_options_tot;
    std::cout << "x: " << x << " y: " << y << " sel: " << select << std::endl;
    if (select < max_ip)
    {
        //pick a bond from in-plane neighbors
        select_move(x, y, select);
    }
    if (select >= max_ip)
    {
        move_op(x, y); //bond to out-of-plane neighbors i.e. fill the vacancy
    }
}

void select_move(int x, int y, int sel)
{
    int ip_loc = 0;

    int xx = 0;
    int yy = 0;

    if (ip_loc == sel)
    {
        xx = x;
        yy = y - 1;
        std::cout << "Moved ( " << x << "," << y << ") to A" << std::endl;
    }
    ip_loc++;

    if (ip_loc == sel)
    {

        xx = x + 1;
        yy = y;
        std::cout << "Moved ( " << x << "," << y << ") to B" << std::endl;
    }
    ip_loc++;

    if (ip_loc == sel)
    {
        xx = x;
        yy = y + 1;
        std::cout << "Moved ( " << x << "," << y << ") to C" << std::endl;
    }
    ip_loc++;

    if (ip_loc == sel)
    {
        xx = x - 1;
        yy = y;

        std::cout << "Moved ( " << x << "," << y << ") to D" << std::endl;
    }

    add_elem(xx, yy);
    rem_elem(x, y);
    std::cout << "new height current elem: " << plant[x][y].height << std::endl;
    std::cout << "new height prev elem: " << plant[xx][yy].height << std::endl;
}

void move_op(int x, int y)
{
    plant[x][y].vacant = false;
    std::cout << "Adatom adsorbed" << std::endl;
}

void add_nn(int flux)
{
    for (int idx_nn = 0; idx_nn < flux; idx_nn++)
    {
        int xx = rand() % x_max;
        int yy = rand() % y_max;
        add_elem(xx, yy);
    }
}

void add_elem(int x, int y)
{
    if (plant[x][y].occupied == false)
    {
        plant[x][y].occupied = true;
    }

    (plant[x][y].height)++;
}

void rem_elem(int x, int y)
{
    if (plant[x][y].occupied == false or plant[x][y].height == 0)
    {
        std::cout << "ERROR, removing occupied element at (" << x << "," << y << ")" << std::endl;
        while (1)
        {
        }
    }

    (plant[x][y].height)--;
    if (plant[x][y].height == 0)
    {
        plant[x][y].occupied = false;
    }
}

void print_substrate()
{
    for (int idx_x = 0; idx_x < x_max; idx_x++)
    {
        for (int idx_y = 0; idx_y < y_max; idx_y++)
        {
            if (plant[idx_x][idx_y].occupied == true)
            {
                std::cout << plant[idx_x][idx_y].height << " ";
            }
            else
            {
                std::cout << " X ";
            }
        }
        std::cout << std::endl;
    }
}

void analyze_surface()
{
    float mean = 0.0f;
    float norm = static_cast<float>(x_max * y_max);

    for (int idx_x = 0; idx_x < x_max; idx_x++)
    {
        for (int idx_y = 0; idx_y < y_max; idx_y++)
        {
            mean += float(plant[idx_x][idx_y].height);
        }
    }
    mean = mean / norm;

    float var = 0.0f;
    for (int idx_x = 0; idx_x < x_max; idx_x++)
    {
        for (int idx_y = 0; idx_y < y_max; idx_y++)
        {
            float h = static_cast<float>(plant[idx_x][idx_y].height);
            float dev = std::abs(h - mean);
            var += std::sqrt(2.0f * dev);
        }
    }
    //normalize by number of elements
    var = var / norm;

    std::cout << "mean height: " << mean << std::endl;
    std::cout << "surface roughness: " << var << std::endl;
}