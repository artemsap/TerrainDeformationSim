using System.Collections;
using System.Collections.Generic;
using System.Diagnostics.Contracts;
using UnityEngine;
using System.IO;
using UnityEngine.Analytics;
using static UnityEditor.Experimental.AssetDatabaseExperimental.AssetDatabaseCounters;
using System.Runtime.InteropServices;

struct node_classification
{
    public bool border;
    public int cnt;
}

public class TerrainDeform : MonoBehaviour
{
    public CustomRenderTexture heightMap;
    public Texture2D output;

    public float Kc = 2370.0f;
    public float Kf = 60300.0f;
    public float n = 0.63f;

    public float[,] pressure;

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
        DeformByIntersection();
        //DeformTerrainByObjects(collision);
    }

    private void OnCollisionStay(Collision collision)
    {
        DeformByIntersection();
        //DeformTerrainByObjects(collision);
    }

    public void DeformByIntersection()
    {
        var heights_ter = terrain.terrainData.GetHeights(0, 0, terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution);
        float[,] delta_heights = new float[terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution];
        node_classification[,] classification = new node_classification[terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution];
        Vector3 start = terrain.GetComponent<Transform>().position;
        Vector3 end = terrain.GetComponent<Transform>().position + new Vector3(terrain.terrainData.size.x, 0, terrain.terrainData.size.z);

        float step_x = (float)((end.x - start.x) / (terrain.terrainData.heightmapResolution - 1));
        float step_z = (float)((end.z - start.z) / (terrain.terrainData.heightmapResolution - 1));

        for (int i = 0; i < terrain.terrainData.heightmapResolution; i++)
        {
            for (int j = 0; j< terrain.terrainData.heightmapResolution; j++)
            {
                delta_heights[j, i] = 0;

                //Vector3 upDist = transform.TransformDirection(Vector3.up) * default_height[j, i] * terrain.terrainData.heightmapScale.y;
                RaycastHit hit;
                if (Physics.Raycast(new Vector3(step_x * i, 0, step_z * j), 
                    transform.TransformDirection(Vector3.up),
                    out hit,
                    heights_ter[j, i] * terrain.terrainData.heightmapScale.y))
                {
                    /*Debug.DrawRay(new Vector3(step_x * i, 0, step_z * j), 
                        transform.TransformDirection(Vector3.up) * hit.distance, 
                        Color.green, 
                        2);*/

                    delta_heights[j, i] = (heights_ter[j, i] * terrain.terrainData.heightmapScale.y - hit.distance) / terrain.terrainData.heightmapScale.y;
                    //heights_ter[j, i] = heights_ter[j, i] - delta_heights[j, i];
                }
            }
        }

        //path of file
        //string path = Application.dataPath + "/classification.txt";
        //if (!File.Exists(path))
        //    File.WriteAllText(path, "Startfile: \n\n");

        int n_contacts = 0;

        //На основе тех границ, которые мы определили, классифицируем каждую вершину
        for (int i = 1; i < terrain.terrainData.heightmapResolution - 1; i++)
        {
            for (int j = 1; j < terrain.terrainData.heightmapResolution - 1; j++)
            {
                float left = delta_heights[j, i - 1];
                float right = delta_heights[j, i + 1];
                float up = delta_heights[j + 1, i];
                float down = delta_heights[j - 1, i];

                if (left == 0 && right == 0 && up == 0 && down == 0)
                {
                    //File.AppendAllText(path, classification[j, i].cnt.ToString());
                    continue;
                }

                n_contacts++;

                if (left != 0 && !classification[j, i - 1].border)
                    classification[j, i].cnt++;
                if (right != 0 && !classification[j, i + 1].border)
                    classification[j, i].cnt++;
                if (up != 0 && !classification[j + 1, i].border)
                    classification[j, i].cnt++;
                if (down != 0 && !classification[j - 1, i].border)
                    classification[j, i].cnt++;

                //File.AppendAllText(path, classification[j, i].cnt.ToString());  
            }
            //File.AppendAllText(path, "\n");
        }

        //Определеяем площадь "footprint" и длину его контура

        int countur = 0;
        int ploshad = 0;
        for (int i = 1; i < terrain.terrainData.heightmapResolution - 1; i++)
        {
            for (int j = 1; j < terrain.terrainData.heightmapResolution - 1; j++)
            {
                if (classification[j, i].cnt == 0)
                    continue;

                if (classification[j, i].border)
                    countur += classification[j, i].cnt;
                else
                    ploshad += classification[j, i].cnt;
            }
        }

        pressure = new float[terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution];
        for (int i = 0; i < terrain.terrainData.heightmapResolution; i++)
        {
            for (int j = 0; j < terrain.terrainData.heightmapResolution; j++)
            {
                pressure[j, i] = (Kc * (countur / (2 * ploshad)) + Kf)*Mathf.Pow(delta_heights[j,i] * terrain.terrainData.heightmapScale.y, n); 
            }
        }

        float summ_distance = 0.0f;
        for (int i = 0; i < terrain.terrainData.heightmapResolution; i++)
        {
            for (int j = 0; j < terrain.terrainData.heightmapResolution; j++)
            {
                if (delta_heights[j, i] == 0)
                    continue;

                summ_distance += calc_distance(new Vector2(j, i), delta_heights);
            }
        }

        //path of file
      //  string path = Application.dataPath + "/centrality_coef.csv";
       // if (!File.Exists(path))
       //     File.WriteAllText(path, "\n");

        float[,] centrality_coef = new float[terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution];

        for (int i = 0; i < terrain.terrainData.heightmapResolution; i++)
        {
            for (int j = 0; j < terrain.terrainData.heightmapResolution; j++)
            {
                if (delta_heights[j, i] == 0)
                {
                    //File.AppendAllText(path, centrality_coef[j, i].ToString() + ';');
                    continue;
                }
                    
                centrality_coef[j, i] = (n_contacts * calc_distance(new Vector2(j, i), delta_heights)) / summ_distance;

                //File.AppendAllText(path, centrality_coef[j, i].ToString() + ';');  
            }
            //File.AppendAllText(path, "\n");
        }


        float[,] final_pressure = new float[terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution];
        for (int i = 0; i < terrain.terrainData.heightmapResolution; i++)
        {
            for (int j = 0; j < terrain.terrainData.heightmapResolution; j++)
            {
                final_pressure[j, i] = pressure[j, i] * centrality_coef[j, i];
            }
        }

        terrain.terrainData.SetHeights(0, 0, heights_ter);
    }
    float calc_distance(Vector2 coords, float[,] delta_heights)
    {
        float dist = 0;

        for (int i = 0; i < terrain.terrainData.heightmapResolution; i++)
        {
            for (int j = 0; j < terrain.terrainData.heightmapResolution; j++)
            {
                if (coords.x == j && coords.y == i || delta_heights[j, i] == 0)
                    continue;

                dist += (coords.x - j) * (coords.x - j) + (coords.y - i) * (coords.y - i);
            }
        }

        return 1/dist;
    }

    node_classification[,] detect_border(node_classification[,] nodes, float[,] delta_heights)
    {
        //Определили границы
        for (int i = 1; i < terrain.terrainData.heightmapResolution - 1; i++)
        {
            for (int j = 1; j < terrain.terrainData.heightmapResolution - 1; j++)
            {
                nodes[j, i].cnt = 0;
                nodes[j, i].border = false;

                float left = delta_heights[j, i - 1];
                float right = delta_heights[j, i + 1];
                float up = delta_heights[j + 1, i];
                float down = delta_heights[j - 1, i];
                if (left == 0 && right == 0 && up == 0 && down == 0)
                {
                    continue;
                }

                else if ((left != 0 || right != 0 || up != 0 || down != 0) && (delta_heights[j, i] == 0))
                {
                    nodes[j, i].border = true;
                    //heights_ter[j, i] +=  0.005f; <---- оталдочная проверка корректности найденных границ
                }
            }
        }

        return nodes;
    }

    void DeformTerrainByObjects(Collision collision)
    {
        List<Vector2> list_height_pix = GetHeightPix_collide(collision);

        //Находит минимум и максимумы по обоим осям, чтобы считать только маленький "кусочек" координат террейна
        //Сделано для того, чтобы не считывать огромную матрицу высот для оптимизации
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
        max_x++;
        max_y++;

        var heights_ter = terrain.terrainData.GetHeights(min_y, min_x, max_y - min_y, max_x - min_x);

        //Потенциально этот цикл можно будет распараллелить используя шейдер
        //Цикл с деформацией террейна по заданой логике
        foreach (var pix in list_height_pix)
        {
            int index_x = (int)(pix.x - min_x);
            int index_y = (int)(pix.y - min_y);

            var x = heights_ter[index_x, index_y] - 0.0001f; // <--- Непосредственно здесь происходит деформация

            var min = default_height[(int)pix.x, (int)pix.y] - 0.001f;
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
