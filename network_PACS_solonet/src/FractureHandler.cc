#include "../include/FractureHandler.h"

FractureHandler::FractureHandler ( const GetPot& dataFile,
                                   const size_type& ID,		// ID numero frattura
                                   const std::string& section ) :
    M_ID(ID), M_data(dataFile, section),
            M_integrationMethodVelocity(M_meshFlat),
            M_integrationMethodPressure(M_meshFlat), M_integrationMethodLinear(
                    M_meshMapped), M_integrationMethodPressureVisualization(
                    M_meshMapped), M_meshFEMVelocity(M_meshFlat),
            M_meshFEMPressure(M_meshFlat), M_meshFEMLinear(M_meshMapped),
            M_meshFEMPressureVisualization(M_meshMapped)
{
    M_levelSet.reset(new LevelSetHandler_Type(dataFile, section));
}

void FractureHandler::init ( )
{
    // Geometric transformation usign primal finite elements type in the fracture
    M_geometricTransformation = bgeot::geometric_trans_descriptor(
            M_data.getMeshType());

    //------------------M_mediumMesh di Gamma sul piano orizzontale -----------------------

    // Build a standard M_mediumMesh, ND

    sizeVector_Type fractureNumberSubdivision(M_data.getSpaceDimension());
    std::fill(fractureNumberSubdivision.begin(),
            fractureNumberSubdivision.end(), M_data.getSpatialDiscretization());

    getfem::regular_unit_mesh(M_meshFlat, fractureNumberSubdivision,
            M_geometricTransformation);

    bgeot::base_matrix fractureTransformMatrix(M_data.getSpaceDimension(),
            M_data.getSpaceDimension());

    scalarVector_Type fractureLength(3);
    fractureLength [ 0 ] = M_data.getLengthAbscissa();
    fractureLength [ 1 ] = M_data.getLengthOrdinate();
    fractureLength [ 2 ] = M_data.getLengthQuota();

    for ( size_type i = 0; i < M_data.getSpaceDimension(); ++i )
    {
        fractureTransformMatrix(i, i)
                = (i < M_data.getSpaceDimension()) ? fractureLength [ i ] : 1.0;
    }

    // scale the unit M_mediumMesh to [M_mediumLengthAbscissa,M_mediumLengthOrdinate]
    M_meshFlat.transformation(fractureTransformMatrix);
    for ( size_type i = 0; i < M_meshFlat.points().size(); ++i )
    {
        //M_meshFlat.points () [ i ] [ 0 ] = M_meshFlat.points () [ i ] [ 0 ] + M_data.getTranslateAbscissa ();
        M_meshFlat.points() [ i ] [ 0 ] = M_data.meshSpacing(
                M_meshFlat.points() [ i ] [ 0 ]);
    }
    //------------------M_mediumMesh di Gamma in 3D------------------------------------------

    //costruisco la M_mediumMesh n.2 quella di servizio, quella della frattura mappata

    sizeVector_Type ind(M_meshFlat.points().size());
    for ( size_type i = 0; i < M_meshFlat.points().size(); ++i )
    {
        bgeot::base_node P(M_data.getSpaceDimension() + 2);
        P [ 0 ] = M_meshFlat.points() [ i ] [ 0 ];
        P [ 1 ] = M_levelSet->getData()->z_map(P);
        ind [ i ] = M_meshMapped.add_point(P);
    }

    for ( size_type i = 0; i < M_meshFlat.convex_index().size(); ++i )
    {
        std::vector<bgeot::size_type> point(3);
        point [ 0 ] = M_meshFlat.ind_points_of_convex(i) [ 0 ];
        point [ 1 ] = M_meshFlat.ind_points_of_convex(i) [ 1 ];
        point [ 2 ] = M_meshFlat.ind_points_of_convex(i) [ 2 ];
        M_meshMapped.add_convex(M_geometricTransformation, point.begin());
    }

    // Finite elements setting

    // Finite element type for the dual variable in the fracture
    getfem::pfem fractureFETypeVelocity = getfem::fem_descriptor(
            M_data.getFEMTypeVelocity());

    // Integration type for the dual variable in the fracture
    getfem::pintegration_method fractureIntegrationTypeVelocity =
            getfem::int_method_descriptor(M_data.getIntegrationTypeVelocity());

    // Integration method for the dual variable in the fracture
    M_integrationMethodVelocity.set_integration_method(
            M_meshFlat.convex_index(), fractureIntegrationTypeVelocity);

    //  Finite element space for the dual variable
    M_meshFEMVelocity.set_qdim(M_data.getSpaceDimension());
    M_meshFEMVelocity.set_finite_element(M_meshFlat.convex_index(),
            fractureFETypeVelocity);

    // Finite element type for the ...
    getfem::pfem fractureFEMTypeLinear = getfem::fem_descriptor(
            M_data.getFEMTypeLinear());

    // Integration type for the ...
    getfem::pintegration_method fractureIntegrationTypeLinea =
            getfem::int_method_descriptor(M_data.getIntegrationTypeVelocity());

    // Integration method for the ...
    M_integrationMethodLinear.set_integration_method(
            M_meshMapped.convex_index(), fractureIntegrationTypeLinea);

    //  Finite element space for the ...
    M_meshFEMLinear.set_qdim(M_data.getSpaceDimension());

    M_meshFEMLinear.set_finite_element(M_meshMapped.convex_index(),
            fractureFEMTypeLinear);

    // P(K-1) discontinuous FE (pressure)
    // Finite element type for the primal variable in the fracture

    getfem::pfem fractureFETypePressure = getfem::fem_descriptor(
            M_data.getFEMTypePressure());

    // Integration type for the primal variable in the fracture
    getfem::pintegration_method fractureIntegrationTypePressure =
            getfem::int_method_descriptor(M_data.getIntegrationTypePressure());

    // Integration method for the primal variable in the fracture
    M_integrationMethodPressure.set_integration_method(
            M_meshFlat.convex_index(), fractureIntegrationTypePressure);

    //  Finite element space for the primal variable in the fracture
    M_meshFEMPressure.set_finite_element(M_meshFlat.convex_index(),
            fractureFETypePressure);

    // Integration method for the primal variable in the fracture for visualization
    M_integrationMethodPressureVisualization.set_integration_method(
            M_meshMapped.convex_index(), fractureIntegrationTypePressure);

    //  Finite element space for the primal variable in the fracture for visualization
    M_meshFEMPressureVisualization.set_finite_element(
            M_meshMapped.convex_index(), fractureFETypePressure);

    // Allocate data

    // Allocate the vector for M_etaNormalInterpolated
    gmm::resize(M_etaNormalInterpolated, M_meshFEMPressure.nb_dof());
    gmm::clear(M_etaNormalInterpolated);

    // Allocate the vector for M_etaTangentialInterpolated
    gmm::resize(M_etaTangentialInterpolated, M_meshFEMPressure.nb_dof());
    gmm::clear(M_etaTangentialInterpolated);

    // Allocate the vector for M_muNormalInterpolated
    gmm::resize(M_muNormalInterpolated, M_meshFEMPressure.nb_dof());
    gmm::clear(M_muNormalInterpolated);

    // Allocate the vector for M_muTangentialInterpolated
    gmm::resize(M_muTangentialInterpolated, M_meshFEMPressure.nb_dof());
    gmm::clear(M_muTangentialInterpolated);

    // Fill the vectors M_etaNormalInterpolated, M_etaTangentialInterpolated, M_muNormalInterpolated and M_muTangentialInterpolated of the fracture
    for ( size_type i = 0; i < M_meshFEMPressure.nb_dof(); ++i )
    {
        M_etaNormalInterpolated [ i ] = M_data.etaNormalDistribution(
                M_meshFEMPressure.point_of_dof(i)) * M_data.getEtaNormal();

        M_etaTangentialInterpolated [ i ] = M_data.etaTangentialDistribution(
                M_meshFEMPressure.point_of_dof(i)) * M_data.getEtaTangential();

        M_muNormalInterpolated [ i ] = M_data.muNormalDistribution(
                M_meshFEMPressure.point_of_dof(i)) * M_data.getMuNormal();

        M_muTangentialInterpolated [ i ] = M_data.muTangentialDistribution(
                M_meshFEMPressure.point_of_dof(i)) * M_data.getMuTangential();
    }

    // Allocate the vector for the M_mediummesh size in the fracture
    gmm::resize(M_inverseMeshSize, M_meshFEMPressure.nb_dof());
    gmm::clear(M_inverseMeshSize);

    //M_levelSet.init();

    M_meshFlat.region ( FractureHandler::FRACTURE_UNCUT * ( M_ID + 1 ) ).add ( M_meshFlat.convex_index() );

}

void FractureHandler::computeInvH ( const BCHandlerPtr_Type& bcHandler )
{
    // Computing h^(-1) on external boudaries of the fracture.
    // Useful for impose the boundary condition with Nitsche penalisation
    const sizeVector_Type& fractureDirichlet =
            bcHandler->getFractureBC(M_ID)->getDirichlet();
    const size_type shiftFracture = fractureDirichlet.size();
    for ( size_type i = 0; i < shiftFracture; i++ )
    {
        for ( getfem::mr_visitor
                vis(M_meshFlat.region(fractureDirichlet [ i ])); !vis.finished(); ++vis )
        {
            // Select the current dof
            size_type dof =
                    M_meshFEMPressure.ind_basic_dof_of_element(vis.cv()) [ 0 ];

            // Estimate the element size
            const scalar_type meshSize = M_meshFlat.convex_radius_estimate(
                    vis.cv());

            // Compute h^(-1)
            M_inverseMeshSize [ dof ] = 1.0 / meshSize;

        }
    }

    //normalVectorAndMap( M_mediumMeshFEMPressure);
}

void FractureHandler::normalVectorAndMap ( const getfem::mesh_fem& mediumMeshFEMPressure )
{
    // Assign normal vector and the map of the fracture
    for ( size_type i = 0; i < M_meshFEMPressure.nb_dof(); ++i )
    {

        const base_node& node = mediumMeshFEMPressure.point_of_basic_dof(i);
        const bgeot::dim_type& dim = M_data.getSpaceDimension();
        scalarVector_Type magnificationMapFactor =
                M_levelSet->getData()->map_jac(node, dim);

        scalarVector_Type fractureNormal = M_levelSet->getData()->normal_map(
                node, dim);

        M_magnificationMapFactor1.push_back(magnificationMapFactor [ 0 ]);
        M_normal1.push_back(fractureNormal [ 0 ]);
        M_normal2.push_back(fractureNormal [ 1 ]);
    }
}

size_type FractureHandler::setMeshLevelSetFracture ( FractureHandler& otherFracture, size_type& globalIndex )
{
    const size_type otherFractureId = otherFracture.getId();
    size_type numIntersect = 0;
    if ( !M_meshLevelSetIntersect[ otherFractureId ].get() )
    {
        M_meshLevelSetIntersect[ otherFractureId ].reset ( new GFMeshLevelSet_Type ( M_meshFlat ) );
        LevelSetHandlerPtr_Type otherLevelSet = otherFracture.getLevelSet();
        M_levelSetIntersect [ otherFractureId ].reset ( new GFLevelSet_Type ( M_meshFlat, 1, false )  );
        M_levelSetIntersect [ otherFractureId ]->reinit();

        const size_type nbDof = M_levelSetIntersect [ otherFractureId ]->get_mesh_fem().nb_basic_dof();
        for ( size_type d = 0; d < nbDof; ++d )
        {
            base_node node = M_levelSetIntersect [ otherFractureId ]->get_mesh_fem().point_of_basic_dof(d);
            base_node mappedNode ( node.size() + 1 );
            for ( size_type i = 0; i < node.size(); ++i )
            {
                mappedNode [ i ] = node[i];
            }
            mappedNode [ node.size() ] = M_levelSet->getData()->z_map ( node );

            M_levelSetIntersect [ otherFractureId ]->values(0)[d] = otherLevelSet->getData()->levelSetFunction ( mappedNode );
        }

        M_meshLevelSetIntersect[ otherFractureId ]->add_level_set ( *M_levelSetIntersect [ otherFractureId ] );
        M_meshLevelSetIntersect[ otherFractureId ]->adapt ();

        size_type i_cv = 0;
        dal::bit_vector bc_cv = M_meshLevelSetIntersect[ otherFractureId ]->linked_mesh().convex_index();

        for ( i_cv << bc_cv; i_cv != size_type(-1); i_cv << bc_cv )
        {
            if (M_meshLevelSetIntersect[ otherFractureId ]->is_convex_cut ( i_cv ) )
            {
                M_meshFlat.region ( FractureHandler::FRACTURE_UNCUT * ( M_ID + 1 ) ).sup ( i_cv );
                M_meshFlat.region ( FractureHandler::FRACTURE_INTERSECT * ( M_ID + 1 ) + otherFractureId + 1 ).add( i_cv );
                M_extendedPressure.push_back ( M_meshFEMPressure.ind_basic_dof_of_element ( i_cv )[0] );
                M_extendedVelocity.push_back ( M_meshFEMVelocity.ind_basic_dof_of_element ( i_cv )[0] );
                M_extendedVelocity.push_back ( M_meshFEMVelocity.ind_basic_dof_of_element ( i_cv )[1] );

                M_fractureIntersectElements [ otherFractureId ].push_back ( i_cv );

                pairSize_Type coppia;
                coppia.first = globalIndex;
                coppia.second = 0;
                M_fractureIntersectElementsGlobalIndex [ otherFractureId ].push_back ( coppia );
                ++globalIndex;
                ++numIntersect;
            }
        }

    }

    return numIntersect;
}

size_type FractureHandler::getNumIntersections () const
{
    size_type total = 0;

    for ( size_type i = 0; i < M_fractureIntersectElements.size(); ++i )
    {
        total += M_fractureIntersectElements[i].size();
    }

    return total;
}
