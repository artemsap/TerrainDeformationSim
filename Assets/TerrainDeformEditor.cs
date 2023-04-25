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
            targ.DeformByIntersection(new ObjectInfo(targ.customVelocity, targ.customSize, targ.customPosition));
        }

        if (GUILayout.Button("Reset terrain"))
        {
            targ.ResetTerrain();
        }

        if(GUILayout.Button("Read default height"))
        {
            targ.ReadDefaultHeight();
        }
    }
}
