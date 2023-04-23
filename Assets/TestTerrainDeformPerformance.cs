using System.Collections;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;
using UnityEngine;

public class TestTerrainDeformPerformance : MonoBehaviour
{
    public Terrain terrain;
    public int delay_thread;

    public int terrain_size;

    private bool needUpdateTerrain;

    float[,] heightMap;

    // Start is called before the first frame update
    void Start()
    {
        needUpdateTerrain = false;
        heightMap = terrain.terrainData.GetHeights(0, 0, terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution);
        terrain_size = terrain.terrainData.heightmapResolution;

        Thread myThread = new Thread(new ThreadStart(DeformTerrain));
        myThread.Start();
    }

    public void DeformTerrain()
    {
        while(true)
        {
            Thread.Sleep(delay_thread);
            for (int i = 0; i < terrain_size; i++)
            {
                for (int j = 0; j < terrain_size; j++)
                {
                    heightMap[j, i] += 0.0001f;
                }
            }
            needUpdateTerrain = true;
        }
    }

    public void Update()
    {
        if (needUpdateTerrain)
        {
            terrain.terrainData.SetHeights(0, 0, heightMap);
            needUpdateTerrain = false;
        }
    }

}
