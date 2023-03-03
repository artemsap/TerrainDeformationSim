using System.Collections;
using System.Collections.Generic;
using System.Diagnostics.Contracts;
using UnityEngine;

public class TerrainDefform : MonoBehaviour
{
    public CustomRenderTexture heightMap;
    public Texture2D output;

    Terrain terrain;

    Vector3 terrain_pos_global;
    float[,] default_height;

    // Start is called before the first frame update
    void Start()
    {
        //CustomRenderTexture.active = heightMap;
        //output = new Texture2D(heightMap.width, heightMap.height, TextureFormat.R16, false);
        terrain = GetComponent<Terrain>();
        terrain_pos_global = terrain.GetComponent<Transform>().position;
        terrain_pos_global.x += terrain.terrainData.size.x / 2;
        terrain_pos_global.z += terrain.terrainData.size.z / 2;
        default_height = terrain.terrainData.GetHeights(0, 0, terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution);
    }

    void OnApplicationQuit()
    {
        terrain.terrainData.SetHeights(0, 0, default_height);
    }    

    private void OnCollisionEnter(Collision collision)
    {
        //DeformTerrainByObjects(collision);
        DeformTerrainByObjects_Opt(collision);
    }

    private void OnCollisionStay(Collision collision)
    {
        //DeformTerrainByObjects(collision);
        DeformTerrainByObjects_Opt(collision);
    }

    void DeformTerrainByObjects(Collision collision)
    {
        List<Vector2> list_height_pix = GetHeightPix_collide(collision);

        var heights_ter = terrain.terrainData.GetHeights(0, 0, terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution);
        foreach (var pix in list_height_pix)
        {
            heights_ter[(int)pix.x, (int)pix.y] = Mathf.Clamp(heights_ter[(int)pix.x, (int)pix.y] - 0.0001f,
                                                              default_height[(int)pix.x, (int)pix.y] - 0.002f,
                                                              default_height[(int)pix.x, (int)pix.y]);
        }
        terrain.terrainData.SetHeights(0, 0, heights_ter);
    }

    void DeformTerrainByObjects_Opt(Collision collision)
    {
        List<Vector2> list_height_pix = GetHeightPix_collide(collision);

        int min_x = terrain.terrainData.heightmapResolution;
        int min_y = terrain.terrainData.heightmapResolution;
        int max_x = 0, max_y = 0;

        foreach (var pix in list_height_pix)
        {
            if (pix.x < min_x) min_x = (int)pix.x;
            if (pix.y < min_y) min_y = (int)pix.y;
            if (pix.x > max_x) max_x = (int)pix.x;
            if (pix.y > max_y) max_y = (int)pix.y;
        }

        max_x += 1;
        max_y += 1;

        var heights_ter = terrain.terrainData.GetHeights(min_y, min_x, max_y - min_y, max_x - min_x);
        foreach (var pix in list_height_pix)
        {
            int index_x = (int)(pix.x - min_x);
            int index_y = (int)(pix.y - min_y);
            var x = heights_ter[index_x, index_y] - 0.0001f;
            var min = default_height[(int)pix.x, (int)pix.y] - 0.01f;
            var max = default_height[(int)pix.x, (int)pix.y];
            heights_ter[index_x, index_y] = Mathf.Clamp(x, min, max);
        }
        terrain.terrainData.SetHeights(min_y, min_x, heights_ter);
    }

    List<Vector2> GetHeightPix_collide(Collision collision)
    {
        List<Vector2> list_height_pix = new List<Vector2>();

        List<ContactPoint> contacts = new List<ContactPoint>();
        var cont = collision.GetContacts(contacts);

        foreach (ContactPoint contact in contacts)
        {
            if (contact.thisCollider != null)
            {
                //координаты на террейне
                var local_x = (contact.point.x - terrain_pos_global.x) + terrain.terrainData.size.x / 2;
                var local_z = (contact.point.z - terrain_pos_global.z) + terrain.terrainData.size.z / 2;

                //процентное соотношение где соприкасновение произошло относительно лок. координат террейна
                var per_x = local_x / terrain.terrainData.size.x;
                var per_z = local_z / terrain.terrainData.size.z;

                //Переводим координаты на террейне в пиксели у карты высот
                Vector2 height_vert = new Vector2((int)(terrain.terrainData.heightmapResolution * per_z),
                                                  (int)(terrain.terrainData.heightmapResolution * per_x));

                if (!list_height_pix.Contains(height_vert))
                {
                    //Debug.DrawRay(contact.point, new Vector3(0, 1, 0), Color.green, 100);
                    list_height_pix.Add(height_vert);
                }
            }
        }
        return list_height_pix;
    }
}
