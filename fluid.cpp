#include "fluid.h"

namespace FL
{
    Fluid::Fluid(const uint32_t width, const uint32_t height)
    : Width(width),
      Height(height),
      Step(1.0 / (width < height ? static_cast<double>(width) : static_cast<double>(height))),
      VelocityX0(width + 2U, std::vector<double>(height + 2U, 0.0)),
      VelocityX1(width + 2U, std::vector<double>(height + 2U, 0.0)),
      VelocityY0(width + 2U, std::vector<double>(height + 2U, 0.0)),
      VelocityY1(width + 2U, std::vector<double>(height + 2U, 0.0)),
      Density0(width + 2U, std::vector<double>(height + 2U, 0.0)),
      Density1(width + 2U, std::vector<double>(height + 2U, 0.0)),
      VelocityForceX(width + 2U, std::vector<double>(height + 2U, 0.0)),
      VelocityForceY(width + 2U, std::vector<double>(height + 2U, 0.0)),
      Dt(DT_DEFAULT),
      Viscosity(VISCOSITY_DEFAULT),
      DiffusionConstant(DIFFUSION_DEFAULT)
    {
    }

    void Fluid::setDt(const double dt)
    {
      Dt = dt;
    }

    void Fluid::setViscosity(const double viscosity)
    {
      Viscosity = viscosity;
    }

    void Fluid::setDiffusionConstant(const double diffusion)
    {
      DiffusionConstant = diffusion;
    }

    void Fluid::addVelocityForce(const uint32_t i, const uint32_t j, const double xForce, const double yForce)
    {
      VelocityForceX[i + 1U][j + 1U] += xForce;
      VelocityForceY[i + 1U][j + 1U] += yForce;
    }

    void Fluid::addParticles(const uint32_t i, const uint32_t j, const double particles)
    {
      Density0[i + 1U][j + 1U] = particles;
    }

    double Fluid::getDensity(const uint32_t i, const uint32_t j)
    {
      return Density0[i + 1U][j + 1U];
    }

    velocity_t Fluid::getVelocity(const uint32_t i, const uint32_t j)
    {
      velocity_t v;
      v[0U] = VelocityX0[i + 1U][j + 1U];
      v[1U] = VelocityY0[i + 1U][j + 1U];
      return v;
    }

    // Forcer les bords
    void Fluid::setEdges2(matrix& x, matrix& y)
    {
      // Remise à zéro sur l'axe du bord Y
      for (uint32_t j = 1U; j <= Height; ++j)
      {
        x[0U][j] = 0.0;
        y[0U][j] = y[1U][j];

        x[Width + 1U][j] = 0.0;
        y[Width + 1U][j] = y[Width][j];
      }

      // Remise à zéro sur l'axe du bord X
      for (uint32_t i = 1U; i <= Width; ++i)
      {
        x[i][0U] = x[i][1U];
        y[i][0U] = 0.0;

        x[i][Height + 1U] = x[i][Height];
        y[i][Height + 1U] = 0.0;
      }

      // Moyennes sur les quatre coins
      x[0U][0U] = (x[1U][0U] + x[0U][1U]) / 2.0;
      y[0U][0U] = (y[1U][0U] + y[0U][1U]) / 2.0;

      x[0U][Height + 1U] = (x[1U][Height + 1U] + x[0U][Height]) / 2.0;
      y[0U][Height + 1U] = (y[1U][Height + 1U] + y[0U][Height]) / 2.0;

      x[Width + 1U][Height + 1U] = (x[Width][Height + 1U] + x[Width + 1U][Height]) / 2.0;
      y[Width + 1U][Height + 1U] = (y[Width][Height + 1U] + y[Width + 1U][Height]) / 2.0;

      x[Width + 1U][0U] = (x[Width][0U] + x[Width + 1U][1U]) / 2.0;
      y[Width + 1U][0U] = (y[Width][0U] + y[Width + 1U][1U]) / 2.0;
    }

    // Forcer les bords
    void Fluid::setEdges1(matrix& m)
    {
      // Recopie valeur sur les bordures Y
      for (uint32_t j = 1U; j <= Height; ++j)
      {
        m[0U][j] = m[1U][j];
        m[Width + 1U][j] = m[Width][j];
      }

      // Recopie valeur sur les bordures X
      for (uint32_t i = 1U; i <= Width; ++i)
      {
        m[i][0U] = m[i][1U];
        m[i][Height + 1U] = m[i][Height];
      }

      // Moyennes sur les quatre coins
      m[0U][0U] = (m[1U][0U] + m[0U][1U]) / 2.0;
      m[0U][Height + 1U] = (m[1U][Height + 1U] + m[0U][Height]) / 2.0;
      m[Width + 1U][Height + 1U] = (m[Width][Height + 1U] + m[Width + 1U][Height]) / 2.0;
      m[Width + 1U][0U] = (m[Width][0U] + m[Width + 1U][1U]) / 2.0;
    }

    // Calcul de la viscosité
    void Fluid::computeViscosity(const matrix& velocityX0, const matrix& velocityY0, matrix& velocityX1, matrix& velocityY1)
    {
      // Taux de viscosité à chaque itération sur la grille
      const double rate = Dt * Viscosity * Width * Height;
      for (int k = 0; k < 20; ++k)
      {
        for (uint32_t i = 1U; i <= Width; ++i)
        {
          for (uint32_t j = 1U; j <= Height; ++j)
          {
            computeComponentRate(rate, i, j, velocityX0, velocityX1);
            computeComponentRate(rate, i, j, velocityY0, velocityY1);
          }
        }
        setEdges2(velocityX1, velocityY1);
      }
    }

    // Calcul de la viscosité par composante 
    void Fluid::computeComponentRate(const double rate, const size_t i, const size_t j, const matrix& m0, matrix& m1)
    {
      // Integration
      m1[i][j] = (m0[i][j] + rate * (m1[i - 1U][j] + m1[i + 1U][j] + m1[i][j - 1U] + m1[i][j + 1U])) / (1 + rate * 4);
    }

    // Calcul de la force
    void Fluid::computeForces(matrix& x, matrix& y)
    {
      for (uint32_t i = 1U; i <= Width; ++i)
      {
        for (uint32_t j = 1U; j <= Height; ++j)
        {
          x[i][j] = x[i][j] + VelocityForceX[i][j] * Dt;
          y[i][j] = y[i][j] + VelocityForceY[i][j] * Dt;
          VelocityForceX[i][j] = 0.0;
          VelocityForceY[i][j] = 0.0;
        }
      }
      setEdges2(x, y);
    }

    // Calcul des coefficients de l'advection
    void Fluid::computeAdvectionCoefficients(const uint32_t i, const uint32_t j, const matrix& x0, const matrix& y0, const double ratio, uint32_t& id, uint32_t& jd, double& aSW, double& aSE, double& aNW, double& aNE)
    {
      double xd = i - ratio * x0[i][j];
      double yd = j - ratio * y0[i][j];
      // Rebascule dans les bords
      if (xd < 0.5) xd = 0.5;
      if (xd > (Width + 0.5)) xd = Width + 0.5;
      if (yd < 0.5) yd = 0.5;
      if (yd > (Height + 0.5)) yd = Height + 0.5;
      // Calcul de id et jd
      id = static_cast<uint32_t>(xd);
      jd = static_cast<uint32_t>(yd);
      // Calcul des aires sur les cellules en chevauchement
      aSW = (1.0 + static_cast<double>(id) - xd) * (1.0 + static_cast<double>(jd) - yd);
      aSE = (xd - static_cast<double>(id)) * (1.0 + static_cast<double>(jd) - yd);
      aNW = (1.0 + static_cast<double>(id) - xd) * (yd - static_cast<double>(jd));
      aNE = (xd - static_cast<double>(id)) * (yd - static_cast<double>(jd));
    }

    // Calcul de l'advection
    void Fluid::computeVelocityAdvection(const matrix& x0, const matrix& y0, matrix& x1, matrix& y1)
    {
      const double ratio = Dt / Step;
      for (uint32_t i = 1U; i <= Width; ++i)
      {
        for (uint32_t j = 1U; j <= Height; ++j)
        {
          // Calcul des coefficients de l'advection
          uint32_t id = 0U; 
          uint32_t jd = 0U;
          double aSW = 0.0;
          double aSE = 0.0; 
          double aNW = 0.0;
          double aNE = 0.0;
          computeAdvectionCoefficients(i, j, x0, y0, ratio, id, jd, aSW, aSE, aNW, aNE);
          // Calcul de la valeur proportionnelle à l'aire de chaque cellule
          x1[i][j] = aSW * x0[id][jd] + aSE * x0[id + 1][jd] + aNW * x0[id][jd + 1] + aNE * x0[id + 1][jd + 1];
          y1[i][j] = aSW * y0[id][jd] + aSE * y0[id + 1][jd] + aNW * y0[id][jd + 1] + aNE * y0[id + 1][jd + 1];
        }
      }
      setEdges2(x1, y1);
    }

    void Fluid::computeDivergenceRate(const matrix& x, const matrix& y, matrix& divRate)
    {
      // Calcul hors bord
      for (uint32_t i = 1U; i <= Width; ++i)
      {
        for (uint32_t j = 1U; j <= Height; ++j)
        {
          divRate[i][j] = (((x[i + 1U][j] - x[i - 1U][j]) + (y[i][j + 1U] - y[i][j - 1U])) / 2.0) * Step;
        }
      }
      // Calcul sur les bords
      for (uint32_t i = 1U; i <= Width; ++i)
      {
        divRate[i][0U] = divRate[i][1U];
        divRate[i][Height + 1U] = divRate[i][Height];
      }
      for (uint32_t j = 1U; j <= Height; ++j)
      {
        divRate[0U][j] = divRate[1U][j];
        divRate[Width + 1U][j] = divRate[Width][j];
      }
      divRate[0U][0U] = (divRate[1U][0U] + divRate[0U][1U]) / 2.0;
      divRate[Width + 1U][0U] = (divRate[Width][0U] + divRate[Width + 1U][1U]) / 2.0;
      divRate[0U][Height + 1U] = (divRate[0U][Height] + divRate[1U][Height + 1U]) / 2.0;
      divRate[Width + 1U][Height + 1U] = (divRate[Width + 1U][Height] + divRate[Width][Height + 1U]) / 2.0;
    }

    // Calcul de la projection
    void Fluid::computeProjection(matrix& x, matrix& y)
    {
      matrix divRate(Width + 2U, std::vector<double>(Height + 2U, 0.0));
      computeDivergenceRate(x, y, divRate);
      
      matrix phi(Width + 2U, std::vector<double>(Height + 2U, 0.0));
      for (int k = 0; k < 20; ++k)
      {
        for (uint32_t i = 1U; i <= Width; ++i)
        {
          for (uint32_t j = 1U; j <= Height; ++j)
          {
            phi[i][j] = (-divRate[i][j] + phi[i - 1U][j] + phi[i + 1U][j] + phi[i][j - 1U] + phi[i][j + 1U]) / 4.0;
          }
        }
        setEdges1(phi);
      }
      for (uint32_t i = 1U; i <= Width; ++i)
      {
        for (uint32_t j = 1U; j <= Height; ++j)
        {
          x[i][j] = x[i][j] - ((phi[i + 1U][j] - phi[i - 1U][j]) / 2.0) / Step;
          y[i][j] = y[i][j] - ((phi[i][j + 1U] - phi[i][j - 1U]) / 2.0) / Step;
        }
      }
      setEdges2(x, y);
    }

    // Calcul de la vitesse
    void Fluid::computeVelocity()
    {
      computeViscosity(VelocityX0, VelocityY0, VelocityX1, VelocityY1);
      computeVelocityAdvection(VelocityX1, VelocityY1, VelocityX0, VelocityY0);
      computeForces(VelocityX0, VelocityY0);
      computeProjection(VelocityX0, VelocityY0);
    }

    // Calcul de l'advection
    void Fluid::computeDensityAdvection(const matrix& x0, const matrix& y0, const matrix& density0, matrix& density1)
    {
      const double ratio = Dt / Step;
      for (uint32_t i = 1U; i <= Width; ++i)
      {
        for (uint32_t j = 1U; j <= Height; ++j)
        {
          // Calcul des coefficients de l'advection
          uint32_t id = 0U; 
          uint32_t jd = 0U;
          double aSW = 0.0;
          double aSE = 0.0; 
          double aNW = 0.0;
          double aNE = 0.0;
          computeAdvectionCoefficients(i, j, x0, y0, ratio, id, jd, aSW, aSE, aNW, aNE);
          // Calcul de la valeur proportionnelle à l'aire de chaque cellule
          density1[i][j] = aSW * density0[id][jd] + aSE * density0[id + 1][jd] + aNW * density0[id][jd + 1] + aNE * density0[id + 1][jd + 1];
        }
      }
      setEdges1(density1);
    } 

    // Calcul de la diffusion
    void Fluid::computeDiffusion(const matrix& density0, matrix& density1)
    {
      const double rate = Dt * DiffusionConstant * Width * Height;
      for (int k = 0; k < 20; ++k)
      {
        for (uint32_t i = 1U; i <= Width; ++i)
        {
          for (uint32_t j = 1U; j <= Height; ++j)
          {
            computeComponentRate(rate, i, j, density0, density1);
          }
        }
        setEdges1(density1);
      }
    }

    // Calcul de la densité
    void Fluid::computeDensity()
    {
      computeDiffusion(Density0, Density1);
      computeDensityAdvection(VelocityX0, VelocityY0, Density1, Density0);
    }
}