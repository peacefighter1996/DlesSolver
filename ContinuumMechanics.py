from ctypes.wintypes import UINT
import sympy as sp
from enum import Enum


class Axis(Enum):
    X = 0
    Y = 1


class Node2d:
    def __init__(self, id, name, x, y):
        self.id = id
        self.name = name
        self.x = x
        self.y = y
        self.ResultForce = sp.zeros(2, 1)

    def ClearForce(self):
        self.ResultForce = sp.zeros(2, 1)


class NodeValue:
    def __init__(self, node: Node2d, axis: Axis, value):
        self.Node = node
        self.Axis = axis
        self.Value = value


class Edge:
    def __init__(self, id: UINT, name: str, nodeStart: Node2d, nodeEnd: Node2d, S):
        self.Id = id
        self.Name = name
        self.NodeStart = nodeStart
        self.NodeEnd = nodeEnd
        self.direction = sp.atan2(
            (nodeEnd.y-nodeStart.y), (nodeEnd.x-nodeStart.x))
        self.S = S
        self.D = sp.Matrix([[-sp.cos(self.direction), -sp.sin(self.direction),
                           sp.cos(self.direction), sp.sin(self.direction)]])
        self.Length = sp.sqrt((nodeEnd.x-nodeStart.x) **
                              2 + (nodeEnd.y-nodeStart.y)**2)
        self.K = sp.Matrix(self.D.transpose() * S/self.Length * self.D)
        self.Stretch = 0
        self.Stress = 0

    def CalculateStretch(self, du):
        deltaU = sp.Matrix([[
            du[(self.NodeStart.id)*2],
            du[(self.NodeStart.id)*2+1],
            du[(self.NodeEnd.id)*2],
            du[(self.NodeEnd.id)*2+1]
        ]])
        self.Stretch = (self.D * deltaU.transpose())[0]
        self.Stress = self.Stretch*self.S/self.Length

        self.ForceVector = self.Stress * \
            sp.Matrix([[sp.cos(self.direction)], [sp.sin(self.direction)]])
        self.NodeStart.ResultForce -= self.ForceVector
        self.NodeEnd.ResultForce += self.ForceVector

        return self.Stretch

    
    

class DlesSolver2D:
    """
    2D Discrete Linear Elastic System Solver
    """

    def __init__(self, nodes: list[Node2d], beams: list[Edge]):
        self.Nodes: list[Node2d] = nodes
        self.SystemMatrix: sp.Matrix = sp.zeros(len(nodes)*2, len(nodes)*2)
        self.Beams: list[Edge] = beams
        for beam in beams:
            self.AddBeam(beam)
        self.Kcc_inv: sp.Matrix = sp.zeros(0, 0)

    def AddBeam(self, beam: Edge):
        nodeStartX = (beam.NodeStart.id)*2
        nodeStartY = (beam.NodeStart.id)*2+1
        nodeEndX = (beam.NodeEnd.id)*2
        nodeEndY = (beam.NodeEnd.id)*2+1

        self.SystemMatrix[nodeStartX, nodeStartX] += beam.K[0]
        self.SystemMatrix[nodeStartX, nodeStartY] += beam.K[1]
        self.SystemMatrix[nodeStartX, nodeEndX] += beam.K[2]
        self.SystemMatrix[nodeStartX, nodeEndY] += beam.K[3]
        self.SystemMatrix[nodeStartY, nodeStartX] += beam.K[4]
        self.SystemMatrix[nodeStartY, nodeStartY] += beam.K[5]
        self.SystemMatrix[nodeStartY, nodeEndX] += beam.K[6]
        self.SystemMatrix[nodeStartY, nodeEndY] += beam.K[7]

        self.SystemMatrix[nodeEndX, nodeStartX] += beam.K[8]
        self.SystemMatrix[nodeEndX, nodeStartY] += beam.K[9]
        self.SystemMatrix[nodeEndX, nodeEndX] += beam.K[10]
        self.SystemMatrix[nodeEndX, nodeEndY] += beam.K[11]
        self.SystemMatrix[nodeEndY, nodeStartX] += beam.K[12]
        self.SystemMatrix[nodeEndY, nodeStartY] += beam.K[13]
        self.SystemMatrix[nodeEndY, nodeEndX] += beam.K[14]
        self.SystemMatrix[nodeEndY, nodeEndY] += beam.K[15]

    def SplitSystemMatrix(self, knownDisplacements: list[NodeValue], knownForces: list[NodeValue], invert: bool = True) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix]:
        self.Koo: sp.Matrix = sp.zeros(
            len(knownDisplacements), len(knownDisplacements))
        self.Kcc: sp.Matrix = sp.zeros(
            self.SystemMatrix.cols-len(knownDisplacements), self.SystemMatrix.cols-len(knownDisplacements))
        self.Koc: sp.Matrix = sp.zeros(
            len(knownDisplacements), self.SystemMatrix.cols-len(knownDisplacements))
        self.Kco: sp.Matrix = sp.zeros(
            self.SystemMatrix.rows-len(knownDisplacements), len(knownDisplacements))

        # Sort and prepare the known matrixes
        knownDisplacements.sort(key=lambda x: x.Node.id * 2 + x.Axis.value)
        knownForces.sort(key=lambda x: x.Node.id * 2 + x.Axis.value)
        self.d0 = sp.zeros(len(knownDisplacements), 1)
        self.fc = sp.zeros(len(knownForces), 1)
        for i in range(len(knownDisplacements)):
            self.d0[i] = knownDisplacements[i].Value
        for i in range(len(knownForces)):
            self.fc[i] = knownForces[i].Value

        self.KnownDisplacements = knownDisplacements
        self.KnownForces = knownForces

        # Back fill correction
        self.id0 = sp.zeros(len(knownDisplacements), 1)
        self.idc = sp.zeros(len(knownForces), 1)
        # Fill item by item the different displacements
        cRowIndex = 0
        oRowIndex = 0
        for i in range(self.SystemMatrix.rows):
            if DlesSolver2D.InKnownDisplacements(i, knownDisplacements):
                KooColIndex = 0
                KocColIndex = 0
                for j in range(self.SystemMatrix.cols):
                    if DlesSolver2D.InKnownDisplacements(j, knownDisplacements):
                        self.Koo[oRowIndex,
                                 KooColIndex] = self.SystemMatrix[i, j]
                        KooColIndex += 1
                    else:
                        self.Koc[oRowIndex,
                                 KocColIndex] = self.SystemMatrix[i, j]
                        KocColIndex += 1
                self.id0[oRowIndex] = i
                oRowIndex += 1

            else:
                KccColIndex = 0
                KcoColIndex = 0
                for j in range(self.SystemMatrix.cols):
                    if DlesSolver2D.InKnownDisplacements(j, knownDisplacements):
                        self.Kco[cRowIndex,
                                 KcoColIndex] = self.SystemMatrix[i, j]
                        KcoColIndex += 1
                    else:
                        self.Kcc[cRowIndex,
                                 KccColIndex] = self.SystemMatrix[i, j]
                        KccColIndex += 1
                self.idc[cRowIndex] = i
                cRowIndex += 1
        if (invert):
            self.Kcc_inv = self.Kcc.inv()
        return self.Koo, self.Kcc, self.Koc, self.Kco

    def Solve(self, fa: sp.Matrix) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
        self.dc = self.Kcc_inv * (self.fc - self.Kco*self.d0)
        self.f0 = self.Koc * self.dc + self.Koo * self.d0
        self.fa = fa
        self.fr = self.f0 - fa

        # build force vector and displacement vectors by back fill correction
        self.f = sp.zeros(self.SystemMatrix.rows, 1)
        self.d = sp.zeros(self.SystemMatrix.rows, 1)
        for i in range(len(self.id0)):
            self.f[self.id0[i]] = self.f0[i]
            self.d[self.id0[i]] = self.d0[i]
        for i in range(len(self.idc)):
            self.f[self.idc[i]] = self.fc[i]
            self.d[self.idc[i]] = self.dc[i]

        for node in self.Nodes:
            node.ClearForce()

        for beam in self.Beams:
            beam.CalculateStretch(self.d)

        return self.dc, self.f0, self.fr

    def InKnownDisplacements(i, Displacements):
        for k in range(len(Displacements)):
            if i == Displacements[k].Node.id*2+Displacements[k].Axis.value:
                return True
        return False
