'''
@ The script is developed by Qingbin Liu under the supervision of Jie Liu
E-mail to Qingbin Liu: liuqingb@mail2.sysu.edu.cn; liuqingb@pku.edu.cn
       to discuss technical details.
--------------------------------------------------------------------------------
The script is to convert the data from Abaqus output database format 
to vtk file format for parallel visualization
Python script performs the following three major steps: 
    1) reading Abaqus output file according to the architecture of ODB; 
    2) data decomposition for parallel visualization; 
    3) writing VTK for Paraview.
'''

#import necessary modules to handle Abaqus output database, files and string

from odbAccess import *
from textRepr import *
from string import *
from time import *
import sys, os     
import numpy as np

#Function 1
def ConvertOdb2Vtk(filename = './Odb2Vtk.txt'):  #Modify the default value of filename here to specify the default configuration file
    
    starttime = time()
    # Check if filename points to existing file
    if not os.path.isfile(filename): 
        print 'Parameter file "%s" not found'%(filename)
        sys.exit(2) 
    
    #read odb2vtk file to get parameter setting 
    odb2vtk = open(filename,'rt')
    read = odb2vtk.read()
    input = read.split("'")
    #get odb file's path
    odb_path = input[1]
    #get odb file name
    odbname = input[3]
    #get the output files' path
    vtk_path = input[5]
    #get the quantity of pieces to partition
    numBlocks = int(input[7])
    #get the frame
    outputFrameIndices = (input[9].split("-"))

    input_SDV = input[11].split(", ")

    #get the step
    outputStepIndices = input[13].split(", ")
    #get the instance
    outputInstanceIndices = input[15].split(",")

    #end reding and close odb2vtk file
    odb2vtk.close()
    #display the reading result of odb2vtk file
    print 'odb2vtk reading finished, time elapsed: ', time()-starttime
    print 'Basic Information:'
    print 'Model:',odbname
    print '1 vtk split to', numBlocks, 'blocks'
    print 'Output frames:', outputFrameIndices
    print 'Output Steps:', outputStepIndices
    print 'Output Instances:', outputInstanceIndices
    
    #open an ODB ( Abaqus output database )
    odb = openOdb(os.path.join(odb_path,odbname)+'.odb',readOnly=True)
    print "ODB opened"

    #access geometry and topology information ( odb->rootAssembly->instances->(nodes, elements) )
    rootassembly = odb.rootAssembly
    instances = rootassembly.instances
    #access attribute information
    steps = odb.steps
    #get instance & step information : Quantity and all names
    keyInstances = instances.keys()    
    numInstances = len(keyInstances)
    
    keySteps = odb.steps.keys()
    numSteps = len(keySteps)
    
    # Loop over instances
    for instanceIndex in outputInstanceIndices:
        try:
            instance = instances[instances.keys()[int(instanceIndex)]]
        except IndexError:
            instance = instances[instances.keys()[-1]]
            pass
        print "Instance: ", instance.name
        
        #access nodes & elements
        nodes = instance.nodes
        nodeElementCountList = np.zeros(len(nodes))
        nodeLabeltoIndex = []
        maxNodeLabel = nodes[-1].label
        
        for node in nodes:
            if maxNodeLabel < node.label:
                maxNodeLabel = node.label
                nodeLabeltoIndex.append(-1)
            else:
                nodeLabeltoIndex.append(-1)
        
        elements = instance.elements
        nodeIndex = 0
        for element in elements:
            verticesPerElement = len(element.connectivity)
            for label in element.connectivity:
                if nodeLabeltoIndex < 0:
                    nodeLabeltoIndex[label-1]
                if maxNodeLabel < label:
                    maxNodeLabel = label
                nodeIndex = label - 1
                nodeElementCountList[nodeIndex] += 1 
        
        print 'Nodes not build elements', np.where(nodeElementCountList==0)
        nodeElementCountList[nodeElementCountList==0] = 1
        
        #compute the number of element of each block
        numElementsPerBlock = len(elements)/numBlocks
        numElementsLastBlock = len(elements) - numElementsPerBlock * (numBlocks - 1)
        
    # Loop over steps
    for iStep in outputStepIndices:
    
        try:
            step = odb.steps[odb.steps.keys()[int(iStep)]]
        except IndexError:
            step = odb.steps[odb.steps.keys()[-1]]
            pass
        
        print 'Step: ', step.name
        
        #access attribute(fieldOutputs) information
        frames = step.frames
        
        # Loop over Frames
        for iFrame in outputFrameIndices:
            
            #Detect whether the input frame is out of range
            try:
                currentFrame = frames[int(iFrame)]
            except IndexError:
                currentFrame = frames[-1]
                pass
            
            #Access field output data in each frame
            print 'Frame:',iFrame
            
            if currentFrame.frameValue.is_integer():
                print 'step =',step.name, 'frame =',currentFrame.frameId+1, 'time =', currentFrame.frameValue
                time1 = time()
                fieldOutputData = {}                
                fieldOutputDataType = {}
                for key in currentFrame.fieldOutputs.keys():
                    nodeVariable = currentFrame.fieldOutputs[key]
                    # variable position: NODAL or ELEMENT_NODAL or INTEGRATION_POINT or CENTROID
                    # see scripting user's guide -- chapter 9
                    currentPositions = nodeVariable.locations[0].position
                    if currentPositions == INTEGRATION_POINT:
                        nodeVariable = currentFrame.fieldOutputs[key].getSubset(position=ELEMENT_NODAL)
                    nodeVariableValues = nodeVariable.values
                    
                    # variable type: SCALAR or VECTOR or TENSOR
                    # see scripting user's guide -- chapter 9
                    if nodeVariable.type == SCALAR:
                        nComponents = 1
                    else:
                        nComponents = len(nodeVariable.componentLabels)
                    print key, nodeVariable.name, nodeVariable.type, nodeVariable.componentLabels
                    nodeVariableName = ''
                    '''for string in nodeVariable.description.split(' '):
                        nodeVariableName += string'''
                    
                    fieldOutputData[key] = [[0. for j in range(nComponents)] for i in range(len(nodes))]
                    fieldOutputDataType[key] = nodeVariable.type
                    
                    isNodeVisted = [0]*len(nodes)
                    for value in nodeVariableValues:
                        nodeLabel = value.nodeLabel
                        '''TODO: node index could potentially exceed'''
                        isNodeVisted[nodeLabel-1] += 1
                        if nComponents == 1:
                            fieldOutputData[key][nodeLabel-1][0] += value.data
                        else:
                            for ii, value in enumerate(value.data):
                                fieldOutputData[key][nodeLabel-1][ii] += value
                    
                    # average values on nodes are shared by more than 1 element
                    # only for the case that node values are inteplated from INTEGRATION_POINT
                    if currentPositions == INTEGRATION_POINT:
                        fieldOutputData[key] = (np.array(fieldOutputData[key]).T
                                               /nodeElementCountList).T.tolist()
                                        
                    time1 = time()
                    '''#access Stress components
                    try:
                        Stress = currentFrame.fieldOutputs['S']
                        node_Stress = Stress.getSubset(position=ELEMENT_NODAL)
                        fieldValues = node_Stress.values
                        nComponents = len(Stress.componentLabels)
                        print 'Reading S'
                        for valueX in fieldValues :
                            stressTensors[valueX.nodeLabel-1][0] += 1.
                            for ii, value in enumerate(valueX.data):
                                stressTensors[valueX.nodeLabel-1][ii+1] += value
                            stressTensors[valueX.nodeLabel-1][7] += valueX.mises
                    except (IndexError, KeyError):
                        pass    
                    # can first ave
                    #print "Time elapsed: ", time() - time1, "s"

                    time1 = time()
                    #Logarithmic strain components
                    try:
                        Logarithmic_strain = currentFrame.fieldOutputs['LE']
                        node_Logarithmic_strain = Logarithmic_strain.getSubset(position=ELEMENT_NODAL)
                        fieldValues = node_Logarithmic_strain.values
                        nComponents = len(Logarithmic_strain.componentLabels)
                        print 'Reading LE'
                        for valueX in fieldValues :
                            strainTensors[valueX.nodeLabel-1][0] += 1
                            for ii, value in enumerate(valueX.data):
                                strainTensors[valueX.nodeLabel-1][ii+1] += value
                            strainTensors[valueX.nodeLabel-1][7] += valueX.maxPrincipal
                            strainTensors[valueX.nodeLabel-1][8] += valueX.minPrincipal
                    except (IndexError, KeyError):
                        pass    
                    #print "Time elapsed: ", time() - time1, "s"

                    time1 = time()
                    #Plastic strain components
                    try:
                        Plastic_strain = currentFrame.fieldOutputs['PE']
                        node_Plastic_strain = Plastic_strain.getSubset(position=ELEMENT_NODAL)
                        fieldValues = node_Plastic_strain.values  
                        nComponents = len(Plastic_strain.componentLabels)  
                        print 'Reading PE'
                        for valueX in fieldValues :
                            pstrainTensors[valueX.nodeLabel-1][0] += 1
                            for ii, value in enumerate(valueX.data):
                                pstrainTensors[valueX.nodeLabel-1][ii+1] += value
                            pstrainTensors[valueX.nodeLabel-1][7] += valueX.maxPrincipal
                            pstrainTensors[valueX.nodeLabel-1][8] += valueX.minPrincipal
                    except (IndexError, KeyError):
                        pass    
                    #print "Time elapsed: ", time() - time1, "s"

                    time1 = time()
                    #Equivalent plastic strain
                    try:
                        Equivalent_plastic_strain = currentFrame.fieldOutputs['PEEQ']
                        node_Equivalent_plastic_strain = Equivalent_plastic_strain.getSubset(position=ELEMENT_NODAL)
                        fieldValues = node_Equivalent_plastic_strain.values
                        print 'Reading PEEQ'
                        for valueX in fieldValues :
                            pstrainEquiv[valueX.nodeLabel-1][0] += 1
                            pstrainEquiv[valueX.nodeLabel-1][1] += valueX.data
                    except (IndexError, KeyError):
                        pass    
                    #print "Time elapsed: ", time() - time1, "s"

                    time1 = time()
                    #Equivalent plastic strain
                    for i, iSDV in enumerate(input_SDV):
                        try:
                            tmpSDV = currentFrame.fieldOutputs['SDV'+iSDV]
                            nodeTmpSDV = tmpSDV.getSubset(position=ELEMENT_NODAL)
                            fieldValues = nodeTmpSDV.values
                            print 'Reading SDV'+iSDV
                            for valueX in fieldValues :
                                if i == 0:
                                    SDV[valueX.nodeLabel-1][i] += 1
                                SDV[valueX.nodeLabel-1][i+1] += valueX.data
                        except (IndexError, KeyError):
                            pass
                    '''
                print "Partitionning model and writing vtk files ......"
                #piece cycle, to partion the model and create each piece for vtk files        
                for iBlock in range(numBlocks):
                    time1 = time()
                    print "frame:",currentFrame.frameId+1,"; block:",iBlock
                    
                    #create and open a VTK(.vtu) files
                    if(numBlocks > 1):
                        outfile = open(os.path.join(vtk_path,odbname)+'_'+step.name+'_'+instance.name+'f%03d'%int(iFrame)+' '+'p'+str(iBlock)+'.vtu','w')
                    if(numBlocks == 1):
                        outfile = open(os.path.join(vtk_path,odbname)+'_'+step.name+'_'+instance.name+'f%03d'%int(iFrame)+'.vtu','w')
                    
                    #<VTKFile>, including the type of mesh, version, and byte_order
                    outfile.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
                    #<UnstructuredGrid>
                    outfile.write('<UnstructuredGrid>'+'\n')
                    #<Piece>, including the number of points and cells
                    if(iBlock == numBlocks-1):
                        outfile.write('<Piece NumberOfPoints="%d" NumberOfCells="%d">\n'%(len(nodes),numElementsPerBlock))
                    else:
                        outfile.write('<Piece NumberOfPoints="%d" NumberOfCells="%d">\n'%(len(nodes),numElementsLastBlock))

                    
                    print "Writing Nodes ......"
                    #<Points> Write nodes position into vtk files
                    outfile.write('<Points>\n')
                    outfile.write('<DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
                    for node in nodes:
                        for idim in range(3):
                            outfile.write(' %11.8e'%node.coordinates[idim])
                        outfile.write('\n  ')        
                    outfile.write('\n</DataArray>\n</Points>\n')
                    #</Points>
                    
                    print "Writing Results data ......"
                    #<PointData> Write results data into vtk files
                    tensorStrs = ''
                    vectorStrs = ''
                    scalarStrs = ''
                    for i, key in enumerate(sorted(fieldOutputData.keys())):
                        if(fieldOutputDataType[key].getText()[:6] == 'TENSOR'):
                            if tensorStrs == '':
                                tensorStrs += key
                            else:
                                tensorStrs += ','+key
                        elif (fieldOutputDataType[key].getText()[:6] == 'VECTOR'):
                            if vectorStrs == '':
                                vectorStrs += key
                            else:
                                vectorStrs += ','+key
                        else:
                            if scalarStrs == '':
                                scalarStrs += key
                            else:
                                scalarStrs += ','+key

                    outfile.write('<PointData Tensors="%s" Vectors="%s" Scalars="%s">\n'%(tensorStrs, vectorStrs,scalarStrs))
                    
                    for i, key in enumerate(sorted(fieldOutputData.keys())):
                        #<DataArray>
                        nComponents = np.array(fieldOutputData[key]).shape[-1]
                        outfile.write('<DataArray type="Float32" Name="%s" NumberOfComponents="%d" format="ascii">\n'%(key, nComponents))
                        for ivalue in fieldOutputData[key]:
                            for iComponent in ivalue:
                                outfile.write('%11.8e '%iComponent)
                            outfile.write('\n')
                        outfile.write('</DataArray>\n')
                        #</DataArray>
                        
                    outfile.write("</PointData>"+'\n')
                    #</PointData>
                    '''
                    #Spatial displacement, <DataArray>
                    outfile.write('<DataArray type="Float32" Name="Spatial_displacement" NumberOfComponents="3" format="ascii">\n')
                    for i in reop_N:
                        k = nodes[storedNodeList[i]].label-1
                        X,Y,Z = nodeVectors[k][0],nodeVectors[k][1],nodeVectors[k][2]
                        outfile.write('%11.8e %11.8e %11.8e\n'%(X,Y,Z))
                    outfile.write('</DataArray>\n')
                    #</DataArray>
                    
                    #Reaction force
                    outfile.write('<DataArray type="Float32" Name="Reaction_force" NumberOfComponents="3" format="ascii">\n')
                    for i in reop_N:
                        k = nodes[storedNodeList[i]].label-1
                        X,Y,Z = nodeVectors[k][9],nodeVectors[k][10],nodeVectors[k][11]
                        outfile.write('%11.8e %11.8e %11.8e\n'%(X,Y,Z))
                    outfile.write('</DataArray>\n')    
                    #</DataArray>
                    
                    #Stress Mises, <DataArray>
                    outfile.write('<DataArray type="Float32" Name="von_Mises" format="ascii">\n')
                    for i in reop_N:
                        k = nodes[storedNodeList[i]].label-1
                        X = stressTensors[k][7]/stressTensors[k][0]
                        outfile.write('%11.8e\n'%X)
                    outfile.write('</DataArray>\n')
                    #</DataArray>
                    
                    #Solution-dependent state variable, <DataArray>
                    for ii, iSDV in enumerate(input_SDV):
                        outfile.write('<DataArray type="Float32" Name="SDV' + iSDV + '" format="ascii">\n')
                        for i in reop_N:
                            k = nodes[storedNodeList[i]].label-1
                            X = SDV[k][ii+1]/SDV[k][0]
                            outfile.write('%11.8e\n'%X)
                        outfile.write('</DataArray>\n')
                        #</DataArray>'''
                        

                    
                    print "Writing Cells ......"
                    #<Cells> Write cells into vtk files
                    outfile.write('<Cells>\n')
                    #Connectivity
                    outfile.write('<DataArray type="Int32" Name="connectivity" format="ascii">\n')
                    
                    for element in elements:
                        for label in element.connectivity:
                            outfile.write('%d '%(label - 1))
                        outfile.write('\n')
                    outfile.write('</DataArray>'+'\n')
                    
                    #Offsets
                    outfile.write('<DataArray type="Int32" Name="offsets" format="ascii">'+'\n')
                    offset = 0
                    for element in elements:
                        offset += len(element.connectivity)
                        outfile.write('%d\n'%(offset))
                    outfile.write('</DataArray>\n')
                    #Type
                    outfile.write('<DataArray type="UInt8" Name="types" format="ascii">\n')
                    for element in elements:
                        elementTypeAbaqus = element.type
                        if ( elementTypeAbaqus== 'C3D8R'):
                            #Hexahedron
                            elementTypeParaview = 12
                        elif (elementTypeAbaqus == 'C3D4'):
                            #Tetrahedron
                            elementTypeParaview = 10
                        elif (elementTypeAbaqus == 'CPE4R'):
                            #Quadrilateral
                            elementTypeParaview = 9
                        elif (elementTypeAbaqus == 'CPE3'):
                            #Triangle
                            elementTypeParaview = 5
                        else:
                            print 'Error: element formate not supported!'
                        outfile.write(str(elementTypeParaview)+'\n')
                    outfile.write('</DataArray>\n')
                    outfile.write('</Cells>\n')
                    #</Cells>

        
                    #</Piece>
                    outfile.write('</Piece>'+'\n')
                    #</UnstructuredGrid>
                    outfile.write('</UnstructuredGrid>'+'\n')
                    #</VTKFile>
                    outfile.write('</VTKFile>'+'\n')
                
                    outfile.close()
                    print "Time elapsed: ", time() - time1, "s" 
                
                
            
    odb.close()
    print "Total time elapsed: ", time() - starttime, "s"

