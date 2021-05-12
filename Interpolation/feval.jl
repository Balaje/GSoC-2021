module evalAtPoint

using PolygonOps
using Gridap
using Gridap.FESpaces
using Gridap.ReferenceFEs
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Fields
using Gridap.CellData
using FillArrays
using InteractiveUtils



"""
    Define the evaluate function for an FEFunction and Point:
        evaluate(uh::SingleFieldFEFunction, P::Point)
    Uses the PolygonOps package to locate the point inside Polygon

    TODO: Other FESpaces and MultifieldFEFunctions

"""
function evaluate(uh::SingleFieldFEFunction, P::Point)
    #Φ=get_cell_shapefuns(get_fe_space(uh))

    # Find the cell dofs
    U=get_cell_dof_values(uh)

    # Convert the reference FE to the physical FE.
    nloc=length(U[1]) # Number of local dofs
    phyM=fill(zeros(nloc,),length(U)) # Physical Space Basis Functions
    marker=fill(0,length(U)) # Marker array to locate point
    C=get_cell_coordinates(get_triangulation(uh)) # Coordinates
    for m=1:length(U)
        # For triangles
        if(nloc==3)
            M=[1 0 0; 0 1 0; 0 0 1]
            coords=[1 1 1; C[m][1][1] C[m][2][1] C[m][3][1]; C[m][1][2] C[m][2][2] C[m][3][2]]
            phyM[m] = (M*inv(coords))*[1; P[1]; P[2]]
        end
        # For quads
        if(nloc==4)
            M=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
            coords=[1 1 1 1;
                C[m][1][1] C[m][2][1] C[m][3][1] C[m][4][1];
                C[m][1][2] C[m][2][2] C[m][3][2] C[m][4][2];
                C[m][1][1]*C[m][1][2] C[m][2][1]*C[m][2][2] C[m][3][1]*C[m][3][2] C[m][4][1]*C[m][4][2]
                    ]
            phyM[m] = (M*inv(coords))*[1; P[1]; P[2]; P[1]*P[2]]
        end
        # Marker to locate triangle
        marker[m] = inpolygon(P,[C[m]; C[m][1]])
    end
    # Find the value
    val = unique(trunc.(abs.(marker).*(dot.(U,phyM)),digits=10))
    return sum(val)
end

# Dispatch the above Evaluate to Gridap.evaluate
Gridap.evaluate(uh::SingleFieldFEFunction, P::Point) = evaluate(uh, P)


##END MODULE
end
