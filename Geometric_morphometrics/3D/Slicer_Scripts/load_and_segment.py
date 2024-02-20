from pathlib import Path

# need to type in location of your scans. Here it would be the 'Raw_Scans' directory

working_dir = (Path.cwd().parent.parent.parent / 'data' / 'Morphology' / '3D' / 'Raw_Scans' / 'batchfinal').resolve()
NHDR_paths = [f for f in working_dir.glob('*') if f.suffix == '.nhdr']

for NHDR_path in NHDR_paths:
    # NHDR path ## Cannot have any spaces in the path!! ##
    #NHDR_path = r"...\FGFR-Branches_GM\Raw_Scans\chick_ctr_2.nhdr"
    # NHDR_path = working_dir + "\\" + NHDR_file

    # load volume
    print('loading sample ' + NHDR_path.stem)
    loadedVolumeNode = slicer.util.loadVolume(str(NHDR_path))

    # lots of this is from https://gist.github.com/lassoan/5ad51c89521d3cd9c5faf65767506b37 and
    #                      https://gist.github.com/lassoan/4d0b94bda52d5b099432e424e03aa2b1

    # Create a new segmentation
    segmentationNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode")
    #segmentationNode.CreateDefaultDisplayNodes() # only needed for display
    segmentationNode.SetReferenceImageGeometryParameterFromVolumeNode(loadedVolumeNode)

    # Create temporary segment editor to get access to effects
    segmentEditorWidget = slicer.qMRMLSegmentEditorWidget()
    segmentEditorWidget.setMRMLScene(slicer.mrmlScene)
    segmentEditorNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentEditorNode")
    segmentEditorWidget.setMRMLSegmentEditorNode(segmentEditorNode)
    segmentEditorWidget.setSegmentationNode(segmentationNode)
    segmentEditorWidget.setMasterVolumeNode(loadedVolumeNode)

    # Create a segment 'AutoThresh' and apply autothresholding on it
    print('Thresholding...')
    addedSegmentID = segmentationNode.GetSegmentation().AddEmptySegment('AutoThresh')
    segmentEditorNode.SetSelectedSegmentID(addedSegmentID)

    segmentEditorWidget.setActiveEffectByName("Threshold")
    effect = segmentEditorWidget.activeEffect()
    effect.self().autoThreshold("MOMENTS", "SET_LOWER_MAX")
    effect.self().onApply()

    # Apply morphological smoothing filters
    segmentEditorWidget.setActiveEffectByName("Smoothing")
    effect = segmentEditorWidget.activeEffect()
    effect.setParameter("SmoothingMethod", "MORPHOLOGICAL_CLOSING")
    effect.setParameter("KernelSizeMm", "3.0")
    effect.self().onApply()

    effect.setParameter("SmoothingMethod", "MEDIAN")
    effect.setParameter("KernelSizeMm", "3.0")
    effect.self().onApply()

    effect.setParameter("SmoothingMethod", "MORPHOLOGICAL_OPENING")
    effect.setParameter("KernelSizeMm", "3.0")
    effect.self().onApply()

    effect.setParameter("SmoothingMethod", "MEDIAN")
    effect.setParameter("KernelSizeMm", "3.0")
    effect.self().onApply()

    # Find largest cavity
    segmentEditorWidget.setActiveEffectByName("Islands")
    effect = segmentEditorWidget.activeEffect()
    effect.setParameterDefault("Operation", "KEEP_LARGEST_ISLAND")
    effect.self().onApply()

    # Delete temporary segment editor
    segmentEditorWidget = None
    slicer.mrmlScene.RemoveNode(segmentEditorNode)

    # Save
    print('Saving...')
    subject_name = Path(NHDR_path).stem
    # Save labelmap
    labelmapVolumeNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLabelMapVolumeNode')
    slicer.modules.segmentations.logic().ExportVisibleSegmentsToLabelmapNode(segmentationNode, labelmapVolumeNode, loadedVolumeNode)
    save_filename = subject_name + '_autothresh_label.nrrd'
    save_name = Path(working_dir / save_filename).resolve()
    slicer.util.saveNode(labelmapVolumeNode, str(save_name))

    # Save segmentation
    save_filename = subject_name + '_autothresh.seg.nrrd'
    save_name = Path(working_dir / save_filename).resolve()
    slicer.util.saveNode(segmentationNode, str(save_name))

    # Save model
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    exportFolderItemId = shNode.CreateFolderItem(shNode.GetSceneItemID(), "Segments")
    slicer.modules.segmentations.logic().ExportAllSegmentsToModels(segmentationNode, exportFolderItemId)
    segmentModels = vtk.vtkCollection()
    shNode.GetDataNodesInBranch(exportFolderItemId, segmentModels)
    modelNode = segmentModels.GetItemAsObject(0)
    save_filename = subject_name + '_autothresh_model.ply'
    save_name = Path(working_dir / save_filename).resolve()
    slicer.util.saveNode(modelNode, str(save_name))

    # Cleanup
    slicer.mrmlScene.RemoveNode(segmentationNode)
    slicer.mrmlScene.RemoveNode(labelmapVolumeNode)
    slicer.mrmlScene.RemoveNode(shNode)
    slicer.mrmlScene.RemoveNode(modelNode)
    slicer.mrmlScene.Clear(0)
print('Done!')