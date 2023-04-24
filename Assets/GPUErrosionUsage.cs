using System.Collections;
using System.Collections.Generic;
using System.Security.Cryptography;
using System.Threading;
using UnityEngine;
using static UnityEditor.PlayerSettings;

public class GPUErrosionUsage : MonoBehaviour
{
    [HideInInspector]
    public ComputeShader ErrosionShader;

    public Terrain terrain;
    public int delay_thread = 100;

    public int terrain_size;

    private bool needUpdateTerrain;

    float[,] heightMap;
    int KerID;

    public int fi = 30;

    ComputeBuffer output;

    float[,] output_heights;

    float dzlim;

    // Start is called before the first frame update
    void Start()
    {
        needUpdateTerrain = false;

    }

    private void Update()
    {
        terrain_size = terrain.terrainData.heightmapResolution;

        float res = (terrain.terrainData.size.x / terrain_size);
        dzlim = (res * Mathf.Tan(Mathf.Deg2Rad * fi)) / terrain.terrainData.heightmapScale.y;

        output = new ComputeBuffer(terrain_size * terrain_size, sizeof(float));
        heightMap = terrain.terrainData.GetHeights(0, 0, terrain_size, terrain_size);
        output.SetData(heightMap);

        KerID = ErrosionShader.FindKernel("Errosion");
        ErrosionShader.SetBuffer(KerID, "output_terrain", output);
        ErrosionShader.SetInt("terrain_size", terrain_size);
        ErrosionShader.SetFloat("dzlim", dzlim);

        ErrosionShader.Dispatch(KerID, terrain_size / 8, terrain_size / 8, 1);
        output.GetData(heightMap);

        terrain.terrainData.SetHeights(0, 0, heightMap);
    }

    public void GpuErrosion()
    {
        terrain_size = terrain.terrainData.heightmapResolution;

        float res = (terrain.terrainData.size.x / terrain_size);
        dzlim = (res * Mathf.Tan(Mathf.Deg2Rad * fi)) / terrain.terrainData.heightmapScale.y;

        output = new ComputeBuffer(terrain_size * terrain_size, sizeof(float));
        heightMap = terrain.terrainData.GetHeights(0, 0, terrain_size, terrain_size);
        output.SetData(heightMap);

        KerID = ErrosionShader.FindKernel("Errosion");
        ErrosionShader.SetBuffer(KerID, "output_terrain", output);
        ErrosionShader.SetInt("terrain_size", terrain_size);
        ErrosionShader.SetFloat("dzlim", dzlim);

        ErrosionShader.Dispatch(KerID, terrain_size / 8, terrain_size / 8, 1);
        output.GetData(heightMap);
       
        terrain.terrainData.SetHeights(0, 0, heightMap);
    }

    public void CpuErrosion()
    {
        terrain_size = terrain.terrainData.heightmapResolution;
        heightMap = terrain.terrainData.GetHeights(0, 0, terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution);
        output_heights = new float[terrain_size, terrain_size];
        output_heights = terrain.terrainData.GetHeights(0, 0, terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution);

        float[,] input_terrain = terrain.terrainData.GetHeights(0, 0, terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution);

        float res = (terrain.terrainData.size.x / terrain_size);

        dzlim = (res * Mathf.Tan(Mathf.Deg2Rad * fi)) / terrain.terrainData.heightmapScale.y;

        for (int i = 1; i < terrain.terrainData.heightmapResolution - 1; i++)
        {
            for (int j = 1; j < terrain.terrainData.heightmapResolution - 1; j++)
            {
                float dz_1 = output_heights[i, j] - output_heights[i + 1, j];
                float dz_2 = output_heights[i, j] - output_heights[i - 1, j];
                float dz_3 = output_heights[i, j] - output_heights[i, j + 1];
                float dz_4 = output_heights[i, j] - output_heights[i, j - 1];

                float sum = 0;
                float max_dz = 0;

                if (dz_1 > 0)
                {
                    sum += dz_1;
                    max_dz = Mathf.Max(max_dz, dz_1);
                }
                    
                if (dz_2 > 0)
                {
                    sum += dz_2;
                    max_dz = Mathf.Max(max_dz, dz_2);
                }
       
                if (dz_3 > 0)
                {
                    sum += dz_3;
                    max_dz = Mathf.Max(max_dz, dz_3);
                }
                    
                if (dz_4 > 0)
                {
                    sum += dz_4;
                    max_dz = Mathf.Max(max_dz, dz_4);
                }
                   
                if (sum == 0)
                    continue;
                   
                if (max_dz > dzlim)
                {
                    float step_1 = (max_dz - dzlim) / 2;

                    output_heights[i, j] -= step_1;
                    if (dz_1 > 0)
                        output_heights[i + 1, j] += step_1 * (dz_1 / sum);
                    if (dz_2 > 0)
                        output_heights[i - 1, j] += step_1 * (dz_2 / sum);
                    if (dz_3 > 0)
                        output_heights[i, j + 1] += step_1 * (dz_3 / sum);
                    if (dz_4 > 0)
                        output_heights[i, j - 1] += step_1 * (dz_4 / sum);
                }                
            }
        }

        terrain.terrainData.SetHeights(0, 0, output_heights);
    }

    public void ResetTerrain()
    {
        terrain.terrainData.SetHeights(0, 0, heightMap);
    }
}
