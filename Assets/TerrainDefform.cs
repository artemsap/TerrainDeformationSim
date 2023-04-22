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
    }

    private void OnCollisionStay(Collision collision)
    {
        GameObject obj = collision.gameObject;
        ObjectInformation objinfo = obj.GetComponent<ObjectInformation>();

        DeformByIntersection(objinfo.ObjectVelocity);
    }

    public void DeformByIntersection(Vector3 velocity)
    {
        if (velocity.x == 0 && velocity.y == 0 && velocity.z == 0)
            return;

        //костыль для корректного определения скорости со столкнувшимся объектом
        float tmp = velocity.x;
        velocity.x = velocity.z;
        velocity.z = tmp;
        velocity = -velocity;

        velocity /= velocity.magnitude;

        var heights_ter = terrain.terrainData.GetHeights(0, 0, terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution);
        float[,] delta_heights_final = new float[terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution];

        List<Vector3> nodes_coordinates = new List<Vector3>();

        float[,] delta_heights = DetectIntersection(ref nodes_coordinates, heights_ter);

        float[,] final_delta_sink = new float[terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution];
        float[,] final_delta_buld = new float[terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution];
        BaseDeform(ref delta_heights_final, ref final_delta_sink, ref final_delta_buld, velocity, delta_heights);

        VolumeDistibution(ref delta_heights_final, delta_heights, final_delta_sink, final_delta_buld, velocity);

        ErosionAlgorithm(ref delta_heights_final, delta_heights, nodes_coordinates);

        coef_correctness = 0.0f;
        for (int i = 1; i < terrain.terrainData.heightmapResolution - 1; i++)
        {
            for (int j = 1; j < terrain.terrainData.heightmapResolution - 1; j++)
            {
                heights_ter[j, i] += delta_heights_final[j, i];
                coef_correctness += (heights_ter[j, i] - default_height[j, i]);
            }
        }

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

    public void ResetTerrain()
    {
        terrain.terrainData.SetHeights(0, 0, default_height);
        coef_correctness = 0;
    }

    float[,] DetectIntersection(ref List<Vector3> nodes_coordinates, float[,] heights_ter)
    {
        float[,] delta_heights = new float[terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution];
        Vector3 start = terrain.GetComponent<Transform>().position;
        Vector3 end = terrain.GetComponent<Transform>().position + new Vector3(terrain.terrainData.size.x, 0, terrain.terrainData.size.z);

        float step_x = (float)((end.x - start.x) / (terrain.terrainData.heightmapResolution - 1));
        float step_z = (float)((end.z - start.z) / (terrain.terrainData.heightmapResolution - 1));

        for (int i = 0; i < terrain.terrainData.heightmapResolution; i++)
        {
            for (int j = 0; j < terrain.terrainData.heightmapResolution; j++)
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
                    if (delta_heights[j, i] == 0)
                        continue;

                    nodes_coordinates.Add(new Vector3(j, delta_heights[j, i], i));
                }
            }
        }
        return delta_heights;
    }

    void BaseDeform(ref float[,] delta_heights_final, ref float[,] final_delta_sink, ref float[,] final_delta_buld, Vector3 Velocity, float[,] delta_heights)
    {
        for (int i = 0; i < terrain.terrainData.heightmapResolution - 1; i++)
        {
            for (int j = 0; j < terrain.terrainData.heightmapResolution - 1; j++)
            {
                if (delta_heights[j, i] == 0)
                    continue;

                float chislitel = Velocity.y * Velocity.y;
                float znamenatel = Velocity.x * Velocity.x +
                    Velocity.y * Velocity.y +
                    Velocity.z * Velocity.z;

                if (znamenatel == 0)
                    continue;

                final_delta_sink[j, i] = delta_heights[j, i] * (chislitel / znamenatel);

                float chislitel_buld = Velocity.x * Velocity.x + Velocity.z * Velocity.z;

                final_delta_buld[j, i] = delta_heights[j, i] * (chislitel_buld / znamenatel);

                delta_heights_final[j, i] = -final_delta_sink[j, i] - final_delta_buld[j, i];
            }
        }
    }

    void VolumeDistibution(ref float[,] delta_heights_final, float[,] delta_heights, float[,] final_delta_sink, float[,] final_delta_buld, Vector3 Velocity)
    {
        Vector3 start = terrain.GetComponent<Transform>().position;
        Vector3 end = terrain.GetComponent<Transform>().position + new Vector3(terrain.terrainData.size.x, 0, terrain.terrainData.size.z);

        float step_x = (float)((end.x - start.x) / (terrain.terrainData.heightmapResolution - 1));
        float step_z = (float)((end.z - start.z) / (terrain.terrainData.heightmapResolution - 1));

        node_classification[,] classification = new node_classification[terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution];
        List<Vector3> borders_coordinates = new List<Vector3>();
        classification = detect_border(classification, delta_heights, ref borders_coordinates);

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
                    dist[k] = (step_z * j - step_x * borders_coordinates[k].x) * (step_z * j - step_x * borders_coordinates[k].x) +
                              (step_x * i - step_x * borders_coordinates[k].z) * (step_x * i - step_x * borders_coordinates[k].z) +
                              (delta_heights[j, i] * terrain.terrainData.heightmapScale.y - borders_coordinates[k].y * terrain.terrainData.heightmapScale.y) *
                              (delta_heights[j, i] * terrain.terrainData.heightmapScale.y - borders_coordinates[k].y * terrain.terrainData.heightmapScale.y);
                    summ_dist += dist[k];

                    Vector3 vec_dist = new Vector3(step_z * j - step_z * borders_coordinates[k].x,
                    delta_heights[j, i] * terrain.terrainData.heightmapScale.y - borders_coordinates[k].y * terrain.terrainData.heightmapScale.y,
                    step_z * i - step_z * borders_coordinates[k].z);

                    double cosa = (vec_dist.x * Velocity.x +
                        vec_dist.y * Velocity.y +
                        vec_dist.z * Velocity.z) /
                        (vec_dist.magnitude * Velocity.magnitude);

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

                    double distr_sink_buld_fin = distr_sink_buld[k] / summ_distr_sink_buld;
                    coef_check_1 += distr_sink_buld_fin;
                    delta_height_border[k] += distr_sink_buld_fin * final_delta_buld[j, i];
                }
            }
        }

        for (int k = 0; k < borders_coordinates.Count(); k++)
        {
            delta_heights_final[(int)borders_coordinates[k].x, (int)borders_coordinates[k].z] = (float)delta_height_border[k];
        }
    }

    void ErosionAlgorithm(ref float[,] delta_heights_final, float[,] delta_heights, List<Vector3> node_coordinates)
    {
        //Здесь я буду делать алгоритм эрозии
        float res = (terrain.terrainData.size.x / terrain.terrainData.heightmapResolution);
        float dzlim = (res * Mathf.Tan(Mathf.Deg2Rad * fi)) / terrain.terrainData.heightmapScale.y;

        float[,] delta_erosion = new float[terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution];
        for (int k = 0; k < num_erosion_iter; k++)
        {
            for (int i = matrix_erosion_size / 2; i < terrain.terrainData.heightmapResolution - matrix_erosion_size / 2; i++)
            {
                for (int j = matrix_erosion_size / 2; j < terrain.terrainData.heightmapResolution - matrix_erosion_size / 2; j++)
                {
                    if (delta_heights_final[j, i] == 0)
                        continue;

                    bool find = FindVector(new Vector3(j, delta_heights[j, i], i), node_coordinates);

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

                            bool find2 = FindVector(new Vector3(index_1, delta_heights[index_1, index_2], index_2), node_coordinates);
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

                            bool find2 = FindVector(new Vector3(index_1, delta_heights[index_1, index_2], index_2), node_coordinates);
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
    }
}
