using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;

[CustomEditor(typeof(TerrainDeform))]
public class TerrainDeformEditor : Editor
{
    public override void OnInspectorGUI()
    {
        TerrainDeform targ = (TerrainDeform)target;

        DrawDefaultInspector();

        if (GUILayout.Button("Detect intersection"))
        {
            targ.DeformByIntersection();
        }
    }
}
