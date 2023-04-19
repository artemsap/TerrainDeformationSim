using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;

[CustomEditor(typeof(TerrainDefform))]
public class TerrainDefformEditor : Editor
{
    public override void OnInspectorGUI()
    {
        TerrainDefform targ = (TerrainDefform)target;

        DrawDefaultInspector();

        if (GUILayout.Button("Detect intersection"))
        {
            targ.DeformByIntersection();
        }

        if (GUILayout.Button("Reset terrain"))
        {
            targ.ResetTerrain();
        }
    }
}
