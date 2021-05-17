module evalAtPoint

using PolygonOps
using SparseArrays

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
    #coord2map.(C)
    #phyM=lazy_map(*,phyM,phyP)

    for m=1:length(U)
        # For triangles
        if(nloc==3)
            coords=[1 1 1; C[m][1][1] C[m][2][1] C[m][3][1]; C[m][1][2] C[m][2][2] C[m][3][2]]
            phyM[m] = (inv(coords))*[1; P[1]; P[2]]
        end
        # For quads
        if(nloc==4)
            coords=[1 1 1 1;
                C[m][1][1] C[m][2][1] C[m][3][1] C[m][4][1];
                C[m][1][2] C[m][2][2] C[m][3][2] C[m][4][2];
                C[m][1][1]*C[m][1][2] C[m][2][1]*C[m][2][2] C[m][3][1]*C[m][3][2] C[m][4][1]*C[m][4][2]
                    ]
            phyM[m] = (inv(coords))*[1; P[1]; P[2]; P[1]*P[2]]
        end
        # Marker to locate triangle
        marker[m] = inpolygon(P,[C[m]; C[m][1]])
    end
    # Find the value
    val = unique(trunc.(abs.(marker).*(dot.(U,phyM)),digits=10))
    return sum(val)
end

"""
[------- OPTIMISED USING LAZY_MAPS ------]
Define the evaluate function for an FEFunction and Point:
        evaluate(uh::SingleFieldFEFunction, P::Point)
    Uses the PolygonOps package to locate the point inside Polygon

    TODO: Other FESpaces and MultifieldFEFunctions

"""

function evaluate_optim(uh::SingleFieldFEFunction, P::Point)
    U=get_cell_dof_values(uh)
    C=reshape(get_cell_coordinates(get_triangulation(uh)),(length(U),)) # Coordinates

    # Define functions

    # inpolygon returns a large array with 0,1,-1: Could it be optimised?
    # This part takes the maximum time and memory.
    indx=findall(x->abs(x)>0, inpolygon.(P,coord2polygon.(C))) #Find the index

    nloc=length(U[1]) # Number of local dofs
    MapArr=lazy_map(inv, lazy_map(coord2map,C))
    PtArr=Fill(point2map(P,nloc),length(U))

    loc2glob=lazy_map(*, MapArr,PtArr)
    Uvals=lazy_map(dot,loc2glob,U)
    return Uvals[indx[1]] #Simply return the first index.
end

coord2polygon(C) = [C; C[1]];

function coord2map(C)
    nloc=size(C)[1]
    if(nloc==3)
        return [1 1 1;
                C[1][1] C[2][1] C[3][1];
                C[1][2] C[2][2] C[3][2]
                ];
    elseif(nloc==4)
        return [1 1 1 1;
                C[1][1] C[2][1] C[3][1] C[4][1];
                C[1][2] C[2][2] C[3][2] C[4][2];
                C[1][1]*C[1][2] C[2][1]*C[2][2] C[3][1]*C[3][2] C[4][1]*C[4][2]
                ] ;
    end
end

function point2map(P,nloc)
    if(nloc==3)
        return [1; P[1]; P[2]];
    elseif(nloc==4)
        return [1; P[1]; P[2]; P[1]*P[2]]
    end
end

# Dispatch the above Evaluate to Gridap.evaluate
Gridap.evaluate(uh::SingleFieldFEFunction, P::Point) = sum(unique(evaluate_optim(uh,P)))


##END MODULE
end
