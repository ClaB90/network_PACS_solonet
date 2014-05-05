/*
 * MeshHandler.h
 *
 *  Created on: Apr 1, 2011
 *      Author: fumagalli
 */

#ifndef MESHHANDLER_H_
#define MESHHANDLER_H_ 1

#include "Core.h"
#include "UsefulFunctions.h"
#include "FracturesSet.h"
#include "BCHandler.h"

class MeshHandler
{
public:

    enum
    {
        UNCUT_REGION = 100
    };

            MeshHandler ( const GetPot& dataFile,
                          const std::string& sectionDomain = "" );
    /** imposto la mesh
     * nel file letto con getpot c'è una variabile che dice se la mesh va importata
     * o deve essere costruita
     */
    void setUpMesh ( );

    /**
     * aggiungo gli elementi necessari alle regioni non tagliate
     * dove sono presenti fratture arricchisco gli elementi finiti
     * tolgo gli elementi relativi alle fratture e li sposto, la matrice
     * avrà tutta la parte relativa alle fratture in fondo
     * la numerazione è fatta in modo automatico
     */
    void setUpRegions ( const FracturesSetPtr_Type& fracture );

    // definizione degli elementi finiti per velocità e pressione
    void setUpFEM ( );

    /**
     * definisco i vettori dei circocentri degli elementi della mesh,
     * dei punti medi dei lati e delle lunghezze dei lati della triangolazione
     */
    void computeMeshMeasures ( );

    /**
     * Just to see what elements are cut by the level set M_levelSet:
     * esporto per paraview gli elementi tagliati, 1 se tagliati 0 se no
     */
    void
    printCuttedElements ( const std::string& vtkFolder = "vtk/",
                          const std::string& fileName = "CuttedElements" ) const;



    inline size_type getSpatialDiscretization ( ) const
    {
        return M_spatialDiscretization;
    }

    inline scalar_type getInclination ( ) const
    {
        return M_inclination;
    }

    inline scalar_type getLengthAbscissa ( ) const
    {
        return M_lengthAbscissa;
    }

    inline scalar_type getLengthOrdinate ( ) const
    {
        return M_lengthOrdinate;
    }

    inline scalar_type getLengthQuota ( ) const
    {
        return M_lengthQuota;
    }

    inline std::string getMeshType ( ) const
    {
        return M_meshType;
    }

    inline const getfem::mesh& getMesh ( ) const
    {
        return M_mesh;
    }

    inline getfem::mesh& getMesh ( )
    {
        return M_mesh;
    }

    inline getfem::mesh_level_set& getMeshLevelSet ()
    {
        return M_meshLevelSet;
    }

    inline bgeot::dim_type getSpaceDimension ( ) const
    {
        return M_spaceDimension;
    }

    inline std::string getIntegrationTypeVelocity ( ) const
    {
        return M_integrationTypeVector;
    }

    inline std::string getIntegrationTypePressure ( ) const
    {
        return M_integrationTypeScalar;
    }

    inline std::string getFEMTypeVelocity ( ) const
    {
        return M_fEMTypeVector;
    }

    inline std::string getFEMTypePressure ( ) const
    {
        return M_fEMTypeScalar;
    }

    inline const getfem::mesh_fem& getMeshFEMCoefficients ( ) const
    {
        return M_meshFEMCoefficients;
    }

    inline const getfem::mesh_fem& getMeshFEMScalar ( ) const
    {
        return M_meshFEMScalar;
    }

    inline const getfem::mesh_fem& getMeshFEMVector ( ) const
    {
        return M_meshFEMVector;
    }

    inline const getfem::mesh_im& getIntegrationMethodVector ( ) const
    {
        return M_integrationMethodVector;
    }

    inline const getfem::mesh_im& getIntegrationMethodScalar ( ) const
    {
        return M_integrationMethodScalar;
    }

    inline const getfem::mesh_region& getRegion ( const size_type& regionFlag ) const
    {
        return M_mesh.region(regionFlag);
    }

    inline getfem::mesh_region& getRegion ( const size_type& regionFlag )
    {
        return M_mesh.region(regionFlag);
    }

    inline const scalarVector_Type& getMeshSize ( ) const
    {
        return M_meshSize;
    }

    inline const scalarVector_Type& getInverseMeshSize ( ) const
    {
        return M_inverseMeshSize;
    }

    inline const scalar_type& getInverseMeshSizeDOF ( const size_type& dof ) const
    {
        return M_inverseMeshSize [ dof ];
    }

    inline const sizeVector_Type& getExtendedDOFScalar ( const size_type& id ) const
    {
        return M_extendedDOFScalar [ id ];
    }

    size_type getCountExtendedDOFScalar ( const scalar_type& id ) const;

    size_type getCountExtendedDOFVector ( const scalar_type& id ) const;

    inline const size_type& getExtendedDOFScalar ( const size_type& id,
                                                   const size_type& dof ) const
    {
        return M_extendedDOFScalar [ id ] [ dof ];
    }

    inline const sizeVector_Type& getExtendedDOFVector ( const size_type& id ) const
    {
        return M_extendedDOFVector [ id ];
    }

    inline const size_type& getExtendedDOFVector ( const size_type& id,
                                                   const size_type& dof ) const
    {
        return M_extendedDOFVector [ id ] [ dof ];
    }

    size_type getCountExtendedIntersectDOFScalar () const;

    const sizeVector_Type& getExtendedIntersectDOFScalar() const
    {
        return M_extendedIntersectDOFScalar;
    }

    size_type getCountExtendedIntersectDOFVector () const;

    const sizeVector_Type& getExtendedIntersectDOFVector() const
    {
        return M_extendedIntersectDOFVector;
    }

    inline const bgeot::pgeometric_trans& getGeometricTransformation ( ) const
    {
        return M_geometricTransformation;
    }

    inline sizeVector_Type& getNonCut ( )
    {
        return M_nonCut;
    }

    inline const scalar_type& getEdgeLength ( const size_type& dof ) const
    {
        return M_edgeLength [ dof ];
    }

    inline const scalar_type& getCircumcentersDistance ( const size_type& dof ) const
    {
        return M_circumcentersDistance [ dof ];
    }

private:

    // questa sistema la cut region creando una lista dei triangoli che non sono veramente
    // tagliati, in pratica quando A1 o A2 sono minori di una certa tolleranza
    void fixCutRegion ( const FractureHandlerPtr_Type& fracture );

    // the M_mediumMesh
    getfem::mesh M_mesh;

    /// getfem::mesh_level_set --->  Keep informations about a mesh crossed by level-sets.
    getfem::mesh_level_set M_meshLevelSet;

    std::string M_meshType;

    bgeot::dim_type M_spaceDimension;

    //elementi finiti
    // Mesh fem for velocity
    std::string M_fEMTypeVector;
    /** getfem::mesh_fem
     * describe the finite element space on which some variables will be described
     * A mesh, with a set of FEM attached to its convexes is called a mesh_fem
     */
    getfem::mesh_fem M_meshFEMVector;
    // mesh_fem for pressure
    std::string M_fEMTypeScalar;
    getfem::mesh_fem M_meshFEMScalar;
    // mesh_fem for coefficients
    getfem::mesh_fem M_meshFEMCoefficients;

    //metodi di integrazione
    // integration method for vector fields
    /**
     * The finite elements method involves evaluation of integrals of these basis functions
     * (or product of basis functions etc.) on convexes (and faces of convexes). In simple
     * cases (polynomial basis functions and linear geometrical transformation), one can
     * evaluate analytically these integrals. In other cases, one has to approximate it using
     * quadrature formulas. Hence, at each convex is attached an integration method along with
     * the FEM. If you have to use an approximate integration method, always choose carefully
     * its order (i.e. highest degree of the polynomials who are exactly integrated with the
     * method): the degree of the FEM, of the polynomial degree of the geometrical transformation,
     * and the nature of the elementary matrix have to be taken into account. If you are unsure
     * about the appropriate degree, always prefer a high order integration method (which will
     * slow down the assembly) to a low order one which will produce a useless linear-system.
     */
    std::string M_integrationTypeVector;

    /** getfem::mesh_im
     * A mesh, with a set of integration methods attached to its convexes is called a mesh_im
     */
    getfem::mesh_im M_integrationMethodVector;
    // integration method for scalar fields
    std::string M_integrationTypeScalar;
    getfem::mesh_im M_integrationMethodScalar;

    // Geometric transformation usign pressure finite elements type
    bgeot::pgeometric_trans M_geometricTransformation;

    // M_mediumMeshSize = local M_mediumMesh size
    scalarVector_Type M_meshSize;
    // M_mediumInverseMeshSize = 1.0 / M_mediumMeshSize;
    scalarVector_Type M_inverseMeshSize;

    std::string M_meshExternal;
    std::string M_meshFolder;

    size_type M_spatialDiscretization;
    scalar_type M_inclination;	// scalar_type è un reale
    scalar_type M_lengthAbscissa;
    scalar_type M_lengthOrdinate;
    scalar_type M_lengthQuota;

    //gradi di libertà estesi (cioè raddoppiati nei triangoli tagliati)

    // The extended dofs for p
    // Base, one sizeVector for each fracture
    sizeVectorContainer_Type M_extendedDOFScalar;
    // Intersect, one sizeVector for each interseciton type
    sizeVector_Type M_extendedIntersectDOFScalar;

    // The extended dofs for u
    // Base, one sizeVector for each fracture
    sizeVectorContainer_Type M_extendedDOFVector;
    // Intersect, one sizeVector for each intersection type
    sizeVector_Type M_extendedIntersectDOFVector;

    // The circumcenters of all the elements
    scalarVector_Type M_circumcentersAbscissa;
    scalarVector_Type M_circumcentersOrdinate;
    scalarVector_Type M_circumcentersDistance;

    // The midpoints for each edge
    scalarVector_Type M_edgeMidpointAbscissa;
    scalarVector_Type M_edgeMidpointOrdinate;

    // Edge length
    scalarVector_Type M_edgeLength;

    // per sistemare il caso degenere (triangolo non veramente tagliato)
    sizeVector_Type M_nonCut;

};

typedef MeshHandler MeshHandler_Type;
typedef boost::shared_ptr<MeshHandler_Type> MeshHandlerPtr_Type;

#endif /* MESHHANDLER_H_ */
