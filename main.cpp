#include "fluid.h"
#include <iostream>
#include <curses.h>
#include <unistd.h>

int main()
{
    initscr();

    uint32_t lines = LINES;
    uint32_t columns = COLS;

   getmaxyx(stdscr, lines, columns);
   //std::cout << "lines " << lines << ", columns " << columns << std::endl;

    FL::Fluid fluid(columns, lines);
    // explosion
    double particles = 4500.0;
    double force = 100.0;
    double gravity = 10;

    fluid.addParticles((columns / 2U) - 2U, lines - 1U, particles);
    fluid.addParticles((columns / 2U) - 2U, lines - 2U, particles);
    fluid.addParticles((columns / 2U) - 1U, lines - 2U, particles);
    fluid.addParticles((columns / 2U), lines - 2U, particles);
    fluid.addParticles((columns / 2U), lines - 1U, particles);

    fluid.addVelocityForce((columns / 2U) - 2U, lines - 1U, -force, 0.0);
    fluid.addVelocityForce((columns / 2U) - 2U, lines - 2U, -force, -force);
    fluid.addVelocityForce((columns / 2U) - 1U, lines - 2U, 0.0, -force);
    fluid.addVelocityForce((columns / 2U), lines - 2U, force, -force);
    fluid.addVelocityForce((columns / 2U), lines - 1U, force, 0.0);

    while (true)
    {
        for (uint32_t i = 0U; i < columns; ++i)
        {
            for (uint32_t j = 0U; j < lines; ++j)
            {
                fluid.addVelocityForce(i, j, 0.0, gravity);

                move(j, i);
                const double density = fluid.getDensity(i, j);
                if (density > 250.0)
                {
                    addch('@');
                }
                else if (density > 200.0)
                {
                    addch('&');
                }
                else if (density > 150.0)
                {
                    addch('#');
                }
                else if (density > 100.0)
                {
                    addch('%');
                }
                else if (density > 75.0)
                {
                    addch('(');
                }
                else if (density > 50.0)
                {
                    addch('/');
                }
                else if (density > 35.0)
                {
                    addch('*');
                }
                else if (density > 20.0)
                {
                    addch(',');
                }
                else if (density > 10.0)
                {
                    addch('.');
                }
                else 
                {
                    addch(' ');
                }
            }
        }
        
        refresh();
        usleep(100000);
        particles /= 1.2;
        force /= 1.2;

        fluid.computeVelocity();
        fluid.computeDensity();
    }

    endwin();
    return 0;
}
