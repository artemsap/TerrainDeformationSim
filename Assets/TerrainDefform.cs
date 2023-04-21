using System.Collections;
using System.Collections.Generic;
using System.Diagnostics.Contracts;
using UnityEngine;
using System.IO;
using UnityEngine.Analytics;
using static UnityEditor.Experimental.AssetDatabaseExperimental.AssetDatabaseCounters;
using System.Runtime.InteropServices;
using System.Runtime.InteropServices.WindowsRuntime;
using static UnityEditor.Progress;
using Unity.Profiling;
using System.Linq;

struct node_classification
{
    public bool border;
    public int cnt;
}

public class TerrainDefform : MonoBehaviour
{
    public CustomRenderTexture heightMap;
    public Texture2D output;

    public bool onlyDeformation;

    public float Kc = 2370.0f;
    public float Kf = 60300.0f;
    public float n = 0.63f;

    public int fi = 45;
    public int num_erosion_iter = 5;
    public int matrix_erosion_size = 3; //Минимум 3

    public Vector3 customVelocity;

    public Terrain terrain;

    public float coef_correctness = 0.0f;

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
        ReadDefaultHeight();
    }

    public void ReadDefaultHeight()
    {
        default_height = terrain.terrainData.GetHeights(0, 0, terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution);
    }

    void OnApplicationQuit()
    {
        terrain.terrainData.SetHeights(0, 0, default_height);
        coef_correctness = 0;
    }    

    private void OnCollisionEnter(Collision collision)
    {
        GameObject obj = collision.gameObject;
        ObjectInformation objinfo = obj.GetComponent<ObjectInformation>();

        DeformByIntersection(objinfo.ObjectVelocity);
        //DeformTerrainByObjects(collision);
    }

    private void OnCollisionStay(Collision collision)
    {
        GameObject obj = collision.gameObject;
        ObjectInformation objinfo = obj.GetComponent<ObjectInformation>();

        DeformByIntersection(objinfo.ObjectVelocity);
        //DeformTerrainByObjects(collision);
    }

    public void DeformByIntersection(Vector3 velocity)
    {
        //костыль для корректного определения скорости со столкнувшимся объектом
        float tmp = velocity.x;
        velocity.x = velocity.z;
        velocity.z = tmp;
        velocity = -velocity;

        velocity /= velocity.magnitude;

        var heights_ter = terrain.terrainData.GetHeights(0, 0, terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution);
        float[,] delta_heights = new float[terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution];
        float[,] delta_heights_final = new float[terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution];
        node_classification[,] classification = new node_classification[terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution];

        List<Vector3> nodes_coordinates = new List<Vector3>();
        List<Vector3> borders_coordinates = new List<Vector3>();

        Vector3 start = terrain.GetComponent<Transform>().position;
        Vector3 end = terrain.GetComponent<Transform>().position + new Vector3(terrain.terrainData.size.x, 0, terrain.terrainData.size.z);

        float step_x = (float)((end.x - start.x) / (terrain.terrainData.heightmapResolution - 1));
        float step_z = (float)((end.z - start.z) / (terrain.terrainData.heightmapResolution - 1));

        for (int i = 0; i < terrain.terrainData.heightmapResolution; i++)
        {
            for (int j = 0; j< terrain.terrainData.heightmapResolution; j++)
            {
                delta_heights[j, i] = 0;
                delta_heights_final[j, i] = 0;

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
                    nodes_coordinates.Add(new Vector3(j, delta_heights[j, i], i));
                    //heights_ter[j, i] = heights_ter[j, i] - delta_heights[j, i];
                }
            }
        }

        classification = detect_border(classification, delta_heights, ref borders_coordinates);
        
        if (!onlyDeformation)
        {
            final_classifiction(classification, delta_heights);

            //Определеяем площадь "footprint" и длину его контура
            int contur = 0;
            int ploshad = 0;
            calc_footprint_params(classification, ref contur, ref ploshad);

            float[,] pressure = calc_pressure(delta_heights, contur, ploshad);

            float[,] centrality_coef = calc_centrality_distr(delta_heights, nodes_coordinates.Count());

            float[,] final_pressure = calc_final_pressure(centrality_coef, pressure);
        }

        float[,] final_delta_sink = new float[terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution];
        float[,] final_delta_buld = new float[terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution];
        for (int i = 0; i < terrain.terrainData.heightmapResolution - 1; i++)
        {
            for (int j = 0; j < terrain.terrainData.heightmapResolution - 1; j++)
            {
                if (delta_heights[j, i] == 0)
                    continue;

                float chislitel = velocity.y * velocity.y;
                float znamenatel = velocity.x * velocity.x +
                    velocity.y * velocity.y +
                    velocity.z * velocity.z;

                final_delta_sink[j, i] = delta_heights[j, i] * (chislitel / znamenatel);

                float chislitel_buld = velocity.x * velocity.x +
                    velocity.z * velocity.z;

                final_delta_buld[j, i] = delta_heights[j, i] * (chislitel_buld / znamenatel);

                delta_heights_final[j, i] = -final_delta_sink[j, i] - final_delta_buld[j, i];
            }
        }

        float correct_detection2 = 0;
        for (int i = 0; i < terrain.terrainData.heightmapResolution - 1; i++)
        {
            for (int j = 0; j < terrain.terrainData.heightmapResolution - 1; j++)
            {
                correct_detection2 += delta_heights_final[j, i];
            }
        }

        print("Correct detection before sinkage distr= " + correct_detection2);

        //Распределеяем вытесненный объем по границам, чтобы потом с помощью эрозии сгладить 
        float[] borders_volume_distr_sink = new float[borders_coordinates.Count()];
        float[] borders_volume_distr_buld = new float[borders_coordinates.Count()];
        double[] delta_height_border = new double[borders_coordinates.Count()];

        for (int i = 0; i < terrain.terrainData.heightmapResolution - 1; i++)
        {
            for (int j = 0; j < terrain.terrainData.heightmapResolution - 1; j++)
            {
                if (delta_heights[j, i] == 0)
                    continue;

                double coef_check = 0;
                double coef_check_1 = 0;

                double[] dist = new double[borders_coordinates.Count()];
                double[] distr_sink_buld = new double[borders_coordinates.Count()];
                double summ_dist = 0;
                double summ_distr_sink_buld = 0;
                for (int k = 0; k < borders_coordinates.Count(); k++)
                {
                    dist[k] = (step_z * j - step_x*borders_coordinates[k].x) * (step_z * j - step_x * borders_coordinates[k].x) +
                              (step_x * i - step_x*borders_coordinates[k].z) * (step_x * i - step_x * borders_coordinates[k].z) +
                              (delta_heights[j, i] * terrain.terrainData.heightmapScale.y - borders_coordinates[k].y * terrain.terrainData.heightmapScale.y) *
                              (delta_heights[j, i] * terrain.terrainData.heightmapScale.y - borders_coordinates[k].y * terrain.terrainData.heightmapScale.y);
                    summ_dist += dist[k];

                    Vector3 vec_dist = new Vector3(step_z * j - step_z * borders_coordinates[k].x,
                    delta_heights[j, i] * terrain.terrainData.heightmapScale.y - borders_coordinates[k].y * terrain.terrainData.heightmapScale.y,
                    step_z * i - step_z * borders_coordinates[k].z);

                    double cosa = (vec_dist.x * velocity.x +
                        vec_dist.y * velocity.y +
                        vec_dist.z * velocity.z) /
                        (vec_dist.magnitude * velocity.magnitude);

                    distr_sink_buld[k] = 0;
                    if (cosa > 0)
                    {
                        distr_sink_buld[k] = cosa;
                    }

                    summ_distr_sink_buld += distr_sink_buld[k];
                }

                for (int k = 0; k < borders_coordinates.Count(); k++)
                {
                    double distr_sink_coef = dist[k] / summ_dist;//1 / dist[k];
                    coef_check += distr_sink_coef;
                    delta_height_border[k] += distr_sink_coef * final_delta_sink[j, i];

                    double distr_sink_buld_fin = distr_sink_buld[k]/ summ_distr_sink_buld;
                    coef_check_1 += distr_sink_buld_fin;
                    delta_height_border[k] += distr_sink_buld_fin * final_delta_buld[j, i];
                }
            }
        }

        for (int k = 0; k < borders_coordinates.Count(); k++)
        {
            delta_heights_final[(int)borders_coordinates[k].x, (int)borders_coordinates[k].z] = (float)delta_height_border[k];
        }

        float correct_detection1 = 0;
        for (int i = 0; i < terrain.terrainData.heightmapResolution - 1; i++)
        {
            for (int j = 0; j < terrain.terrainData.heightmapResolution - 1; j++)
            {
                correct_detection1 += delta_heights_final[j, i];
            }
        }
    
        print("Correct detection after sinkage distr= " + correct_detection1);

        //Здесь я буду делать алгоритм эрозии
        float res = (terrain.terrainData.size.x / terrain.terrainData.heightmapResolution);
        float dzlim = (res * Mathf.Tan(Mathf.Deg2Rad * fi)) / terrain.terrainData.heightmapScale.y;

        float[,] delta_erosion = new float[terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution];
        for (int k = 0; k < num_erosion_iter; k++)
        {
            for (int i = matrix_erosion_size/2; i < terrain.terrainData.heightmapResolution - matrix_erosion_size/2; i++)
            {
                for (int j = matrix_erosion_size/2; j < terrain.terrainData.heightmapResolution - matrix_erosion_size/2; j++)
                {
                    bool find = FindVector(new Vector3(j, delta_heights[j, i], i), nodes_coordinates);

                    if (find)
                        continue;

                    float sum = 0;
                    for (int l = 0; l < matrix_erosion_size; l++)
                    {
                        for (int p = 0; p < matrix_erosion_size; p++)
                        {
                            int index_1 = j - matrix_erosion_size / 2 + l;
                            int index_2 = i - matrix_erosion_size / 2 + p;
                            if (index_1 == j && index_2 == i)
                                continue;

                            bool find2 = FindVector(new Vector3(index_1, delta_heights[index_1, index_2], index_2), nodes_coordinates);
                            if (find2)
                                continue;

                            sum += 1 - delta_heights_final[index_1, index_2];
                        }
                    }

                    if (sum == 0 || (delta_heights_final[j, i] < dzlim && delta_heights_final[j, i] > -dzlim))
                        continue;

                    float step_1 = (delta_heights_final[j, i] - dzlim) / 2;

                    for (int l = 0; l < matrix_erosion_size; l++)
                    {
                        for (int p = 0; p < matrix_erosion_size; p++)
                        {
                            int index_1 = j - matrix_erosion_size / 2 + l;
                            int index_2 = i - matrix_erosion_size / 2 + p;

                            bool find2 = FindVector(new Vector3(index_1, delta_heights[index_1, index_2], index_2), nodes_coordinates);
                            if (find2)
                                continue;

                            if (index_1 == j && index_2 == i)
                                delta_heights_final[index_1, index_2] -= step_1;
                            else
                                delta_heights_final[index_1, index_2] += step_1 * ((1 - delta_heights_final[index_1, index_2]) / sum);

                        }
                    }
                }
            }
        }

        coef_correctness = 0.0f;
        float correct_detection_final = 0;
        for (int i = 1; i < terrain.terrainData.heightmapResolution - 1; i++)
        {
            for (int j = 1; j < terrain.terrainData.heightmapResolution - 1; j++)
            {
                correct_detection_final += delta_heights_final[j, i];
                heights_ter[j, i] += delta_heights_final[j, i];
                coef_correctness += heights_ter[j, i] - default_height[j, i];
            }
        }

        print("Correct detection final= " + correct_detection_final);

        terrain.terrainData.SetHeights(0, 0, heights_ter);
    }

    bool FindVector(Vector3 vectortofind, List<Vector3> nodes_coordinates)
    {
        bool find = false;
        foreach (var node in nodes_coordinates)
        {
            if (vectortofind == node)
            {
                find = true;
                break;
            }
        }

        return find;
    }

    Vector3 calc_vector(Vector3 a, Vector3 b, Vector3 c)
    {
        Vector3 side_1 = b - a;
        Vector3 side_2 = c - a;

        return Vector3.Cross(side_1, side_2);
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

    node_classification[,] detect_border(node_classification[,] nodes, float[,] delta_heights, ref List<Vector3> borders_coordinates)
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
                    borders_coordinates.Add(new Vector3(j, delta_heights[j, i], i));
                    //heights_ter[j, i] +=  0.005f; <---- оталдочная проверка корректности найденных границ
                }
            }
        }

        return nodes;
    }

    node_classification[,] final_classifiction(node_classification[,] nodes, float[,] delta_heights)
    {
        //string path = Application.dataPath + "/classification.txt";
        //if (!File.Exists(path))
        //    File.WriteAllText(path, "Startfile: \n\n");

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
                    //File.AppendAllText(path, nodes[j, i].cnt.ToString());
                    continue;
                }

                if (left != 0 && !nodes[j, i - 1].border)
                    nodes   [j, i].cnt++;
                if (right != 0 && !nodes[j, i + 1].border)
                    nodes[j, i].cnt++;
                if (up != 0 && !nodes[j + 1, i].border)
                    nodes[j, i].cnt++;
                if (down != 0 && !nodes[j - 1, i].border)
                    nodes[j, i].cnt++;

                //File.AppendAllText(path, nodes[j, i].cnt.ToString());  
            }
            //File.AppendAllText(path, "\n");
        }


        return nodes;
    }

    void calc_footprint_params(node_classification[,] nodes, ref int contur, ref int ploshad)
    {
        for (int i = 1; i < terrain.terrainData.heightmapResolution - 1; i++)
        {
            for (int j = 1; j < terrain.terrainData.heightmapResolution - 1; j++)
            {
                if (nodes[j, i].cnt == 0)
                    continue;

                if (nodes[j, i].border)
                    contur += nodes[j, i].cnt;
                else
                    ploshad += nodes[j, i].cnt;
            }
        }
    }

    float[,] calc_pressure(float[,] delta_heights, int contur, int ploshad)
    {
        float [,] pressure = new float[terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution];
        for (int i = 0; i < terrain.terrainData.heightmapResolution; i++)
        {
            for (int j = 0; j < terrain.terrainData.heightmapResolution; j++)
            {
                pressure[j, i] = (Kc * (contur / (2 * ploshad)) + Kf) * Mathf.Pow(delta_heights[j, i] * terrain.terrainData.heightmapScale.y, n);
            }
        }
        return pressure;
    }

    float[,] calc_centrality_distr(float[,] delta_heights, int n_contacts)
    {
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

        return centrality_coef;
    }

    float[,] calc_final_pressure(float[,] centrality_coef, float[,] pressure)
    {
        float[,] final_pressure = new float[terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution];
        for (int i = 0; i < terrain.terrainData.heightmapResolution; i++)
        {
            for (int j = 0; j < terrain.terrainData.heightmapResolution; j++)
            {
                final_pressure[j, i] = pressure[j, i] * centrality_coef[j, i];
            }
        }
        return final_pressure;
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

    public void ResetTerrain()
    {
        terrain.terrainData.SetHeights(0, 0, default_height);
        coef_correctness = 0;
    }
    
}
