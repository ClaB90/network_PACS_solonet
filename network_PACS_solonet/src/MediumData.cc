#include "../include/MediumData.h"

MediumData::MediumData ( const GetPot& dataFile,
                         const std::string& sectionSolver,
                         const std::string& section,
                         const std::string& sectionDomain ) :
            M_section(section),
            M_sectionDomain(M_section + sectionDomain),
            M_sectionSolver(M_section + sectionSolver),
            M_penaltyVector(dataFile(
                    (M_sectionDomain + "penaltyVelocity").data(), 5.)),
            M_penaltyScalar(dataFile(
                    (M_sectionDomain + "penaltyPressure").data(), 1.)),
            // darcy
            M_invK(dataFile((M_sectionSolver + "invK").data(), 1.)),
            M_invKDistribution11(dataFile(
                    (M_sectionSolver + "invKDist11").data(), "1.")),
           M_invKDistribution12(dataFile(
                    (M_sectionSolver + "invKDist12").data(), "1.")),
           M_invKDistribution22(dataFile(
                    (M_sectionSolver + "invKDist22").data(), "1.")),
            M_exact(dataFile((M_sectionSolver + "solution").data(), "x")),
            M_exactInlet(dataFile((M_sectionSolver + "solutionIn").data(),
                    M_exact.data())), M_exactOutlet(dataFile((M_sectionSolver
                    + "solutionOut").data(), M_exact.data())), M_exactFlux(
                    dataFile((M_sectionSolver + "velocity").data(), "1.")),
            M_source(dataFile((M_sectionSolver + "source").data(), "1.")),
            M_exactInitial(dataFile(
                    (M_sectionSolver + "initialCondition").data(), "0."))
{}

// Exact solution, pressure
scalar_type MediumData::exact ( const base_node& x, const scalar_type& t ) const
{
    M_parser.setString(M_exact);

    M_parser.setVariable("x", x [ 0 ]);
    M_parser.setVariable("y", x [ 1 ]);
    M_parser.setVariable("t", t);

    return M_parser.evaluate();
}

// Exact solution, pressure IN (i.e. level set <0)
scalar_type MediumData::exactInlet ( const base_node& x, const scalar_type& t ) const
{
    M_parser.setString(M_exactInlet);

    M_parser.setVariable("x", x [ 0 ]);
    M_parser.setVariable("y", x [ 1 ]);
    M_parser.setVariable("t", t);

    return M_parser.evaluate();
}

// Exact solution, pressure OUT, i.e. level set >0
scalar_type MediumData::exactOutlet ( const base_node& x, const scalar_type& t ) const
{
    M_parser.setString(M_exactOutlet);

    M_parser.setVariable("x", x [ 0 ]);
    M_parser.setVariable("y", x [ 1 ]);
    M_parser.setVariable("t", t);

    return M_parser.evaluate();
}

// Exact solution, velocity (non ho ancora impostato quella corretta)
scalar_type MediumData::exactFlux ( const base_node& x,
                                    const base_node& n,
                                    const scalar_type& t ) const
{
    M_parser.setString(M_exactFlux);

    M_parser.setVariable("x", x [ 0 ]);
    M_parser.setVariable("y", x [ 1 ]);
    M_parser.setVariable("n1", n [ 0 ]);
    M_parser.setVariable("n2", n [ 1 ]);
    M_parser.setVariable("t", t);

    return M_parser.evaluate();
}

// Exact solution, div(Velocity) -- SET = 0 WITH NO MASS SOURCES/SINKS !
scalar_type MediumData::source ( const base_node& x, const scalar_type& t ) const
{
    M_parser.setString(M_source);

    M_parser.setVariable("x", x [ 0 ]);
    M_parser.setVariable("y", x [ 1 ]);
    M_parser.setVariable("t", t);

    return M_parser.evaluate();
}

//questa bella funzione restituisce un'eventuale modulazione del coefficiente di permeabilità
scalar_type MediumData::invKDistribution11 ( const base_node& x ) const
{
    M_parser.setString(M_invKDistribution11);

    M_parser.setVariable("x", x [ 0 ]);
    M_parser.setVariable("y", x [ 1 ]);

    return M_parser.evaluate();
}

//questa bella funzione restituisce un'eventuale modulazione del coefficiente di permeabilità
scalar_type MediumData::invKDistribution12 ( const base_node& x ) const
{
    M_parser.setString(M_invKDistribution12);

    M_parser.setVariable("x", x [ 0 ]);
    M_parser.setVariable("y", x [ 1 ]);

    return M_parser.evaluate();
}

//questa bella funzione restituisce un'eventuale modulazione del coefficiente di permeabilità
scalar_type MediumData::invKDistribution22 ( const base_node& x ) const
{
    M_parser.setString(M_invKDistribution22);

    M_parser.setVariable("x", x [ 0 ]);
    M_parser.setVariable("y", x [ 1 ]);

    return M_parser.evaluate();
}

// Exact solution, pressure
scalar_type MediumData::exactInitial ( const base_node& x ) const
{
    M_parser.setString(M_exactInitial);

    M_parser.setVariable("x", x [ 0 ]);
    M_parser.setVariable("y", x [ 1 ]);

    return M_parser.evaluate();
}
