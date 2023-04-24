using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;

[CustomEditor(typeof(GPUErrosionUsage))]
public class GPUErrosionUsageEditor : Editor
{
    public override void OnInspectorGUI()
    {
        GPUErrosionUsage targ = (GPUErrosionUsage)target;

        DrawDefaultInspector();

        if (GUILayout.Button("Do Errosion"))
        {
            targ.GpuErrosion();
        }

        if (GUILayout.Button("Do CPU Errosion"))
        {
            targ.CpuErrosion();
        }


        if (GUILayout.Button("Reset Terrain"))
        {
            targ.ResetTerrain();
        }
    }
}
