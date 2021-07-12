#include <array>
#include <vector>
#include <cstdint>

// Simulation de l'écoulement d'un fluide - Humbert Florent - 2008
// http://humbert-florent.ftp-developpez.com/algorithmique/fluide/simulation.pdf
// Stable Fluid - Jos Stam - 1999

namespace FL
{
    using velocity_t = std::array<double,2U>;

    class Fluid
    {
    public:
        Fluid(const uint32_t width, const uint32_t height);

        void setDt(const double dt);

        void setViscosity(const double viscosity);

        void setDiffusionConstant(const double constant);

        void addVelocityForce(const uint32_t i, const uint32_t j, const double xForce, const double yForce);

        void addParticles(const uint32_t i, const uint32_t j, const double particles = 100.0);

        // Calcul de la vitesse
        void computeVelocity();

        // Calcul de la densité
        void computeDensity();

        double getDensity(const uint32_t i, const uint32_t j);

        velocity_t getVelocity(const uint32_t i, const uint32_t j);
        
    private:
        // type champ
        using matrix = std::vector<std::vector<double>>;

        // Forcer les bords (champ de vecteur)
        void setEdges2(matrix& x, matrix& y);

        // Forcer les bords (champ de scalaire)
        void setEdges1(matrix& m);

        // Calcul de la viscosité
        void computeViscosity(const matrix& velocityX0, const matrix& velocityY0, matrix& velocityX1, matrix& velocityY1);

        // Calcul de la viscosité par composante 
        void computeComponentRate(const double rate, const size_t i, const size_t j, const matrix& m0, matrix& m1);

        // Calcul de la force
        void computeForces(matrix& x, matrix& y);

        // Calcul des coefficients de l'advection
        void computeAdvectionCoefficients(const uint32_t i, const uint32_t j, const matrix& x0, const matrix& y0, const double ratio, uint32_t& id, uint32_t& jd, double& aSW, double& aSE, double& aNW, double& aNE);

        // Calcul de l'advection
        void computeVelocityAdvection(const matrix& x0, const matrix& y0, matrix& x1, matrix& y1);    

        // Calcul de la divergence
        void computeDivergenceRate(const matrix& x0, const matrix& y0, matrix& divRate);

        // Calcul de la projection
        void computeProjection(matrix& x, matrix& y);

        // Calcul de l'advection
        void computeDensityAdvection(const matrix& x0, const matrix& y0, const matrix& density0, matrix& density1);   

        // Calcul de la diffusion
        void computeDiffusion(const matrix& density0, matrix& density1);

        // Tailles
        uint32_t Width;
        uint32_t Height;
        double Step;

        // Champ de vecteurs
        matrix VelocityX0;
        matrix VelocityX1;
        matrix VelocityY0;
        matrix VelocityY1;

        // Matrice de densité
        matrix Density0;
        matrix Density1;

        // Paramètres
        double Dt;
        double Viscosity;
        double DiffusionConstant;

        // Forces
        matrix VelocityForceX;
        matrix VelocityForceY;

        // Valeurs par défaut
        static constexpr double DT_DEFAULT = 0.02;
        static constexpr double VISCOSITY_DEFAULT = 0.075;
        static constexpr double DIFFUSION_DEFAULT = 0.075;
    };
}