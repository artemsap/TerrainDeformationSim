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
using System.Threading;
using Assembly_CSharp;
using static TMPro.SpriteAssetUtilities.TexturePacker_JsonArray;
using UnityEngine;
struct node_classification
{
    public bool border;
    public int cnt;
}

public struct ObjectInfo
{
    public Vector3 velocity;
    public Vector3 size;
    public Vector3 position;

    public ObjectInfo(Vector3 Velocity, Vector3 Size, Vector3 Position)
    {
        this.velocity = Velocity;
        this.size = Size;
        this.position = Position;
    }
}


public class TerrainDefform : MonoBehaviour
{
    public bool onlyDeformation;

    public float Kc = 2370.0f;
    public float Kf = 60300.0f;
    public float n = 0.63f;

    public float Ma = 0.2f; //в одном метре столбца 

    public int fi = 45;
    public int num_erosion_iter = 5;
    public int matrix_erosion_size = 3; //Минимум 3

    public Vector3 customVelocity;
    public Vector3 customSize;
    public Vector3 customPosition;

    public Terrain terrain;

    public float coef_correctness = 0.0f;

    Vector3 terrain_pos_global;
    float[,] default_height;
    //float[,] heights_ter;
    public int terrain_size;
    public float size_x;
    public float size_z;
    public Vector3 ter_position;
    public float scale_y;

    float[,] heights_ter;
    float[,] delta_heights;
    ObjectInformation objinfo;
    bool needUpdate;


    [HideInInspector]
    public ComputeShader ErrosionShader;

    ComputeBuffer output;
    ComputeBuffer delta;

    bool Conctact;

    // Start is called before the first frame update
    void Start()
    {
        Conctact = false;
        //CustomRenderTexture.active = heightMap;
        //output = new Texture2D(heightMap.width, heightMap.height, TextureFormat.R16, false);
        terrain = GetComponent<Terrain>();
        terrain_pos_global = terrain.GetComponent<Transform>().position;
        terrain_pos_global.x += terrain.terrainData.size.x / 2;
        terrain_pos_global.z += terrain.terrainData.size.z / 2;
        ter_position = terrain.GetComponent<Transform>().position;
        size_x = terrain.terrainData.size.x;
        size_z = terrain.terrainData.size.z;
        scale_y = terrain.terrainData.heightmapScale.y;
        ReadDefaultHeight();
        terrain_size = terrain.terrainData.heightmapResolution;
    }

    public void FixedUpdate()
    {
        if (Conctact)
            return;
        heights_ter = terrain.terrainData.GetHeights(0, 0, terrain_size, terrain_size);
        ErosionAlgorithmGPU(ref heights_ter, delta_heights, 1);
        terrain.terrainData.SetHeights(0, 0, heights_ter);
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
        objinfo = obj.GetComponent<ObjectInformation>();

        DeformByIntersection(new ObjectInfo(objinfo.ObjectVelocity, objinfo.ObjectSize, objinfo.ObjectPosition));
    }

    private void OnCollisionStay(Collision collision)
    {
        Conctact = true;
        GameObject obj = collision.gameObject;
        objinfo = obj.GetComponent<ObjectInformation>();

        DeformByIntersection(new ObjectInfo(objinfo.ObjectVelocity, objinfo.ObjectSize, objinfo.ObjectPosition));
    }

    private void OnCollisionExit(Collision collision)
    {
        Conctact = false;
    }

    public void DeformByIntersection(ObjectInfo obj_info)
    {
        if (obj_info.velocity.x == 0 && obj_info.velocity.y == 0 && obj_info.velocity.z == 0 || obj_info.velocity.y > 0)
            return;

        //костыль для корректного определения скорости со столкнувшимся объектом
        float tmp = obj_info.velocity.x;
        obj_info.velocity.x = obj_info.velocity.z;
        obj_info.velocity.z = tmp;
        obj_info.velocity = -obj_info.velocity;

        obj_info.velocity /= obj_info.velocity.magnitude;

        //terrain = GetComponent<Terrain>();
        //ter_position = terrain.GetComponent<Transform>().position;
        //terrain_pos_global = ter_position;
        //terrain_pos_global.x += terrain.terrainData.size.x / 2;
        //terrain_pos_global.z += terrain.terrainData.size.z / 2;
        
        //size_x = terrain.terrainData.size.x;
        //size_z = terrain.terrainData.size.z;
        //scale_y = terrain.terrainData.heightmapScale.y;
        //terrain_size = terrain.terrainData.heightmapResolution;

        heights_ter = terrain.terrainData.GetHeights(0, 0, terrain_size, terrain_size);
        float[,] delta_heights_final = new float[terrain_size, terrain_size];

        delta_heights = DetectIntersection(out var nodes_coordinates, heights_ter, obj_info);
        needUpdate = true;
        
        BaseDeform(ref delta_heights_final, out var final_delta_sink, out var final_delta_buld, obj_info.velocity, in delta_heights);

        TerrainAirCorrection(ref final_delta_sink, in heights_ter);

        VolumeDistibution(ref delta_heights_final, in delta_heights, in final_delta_sink, in final_delta_buld, obj_info.velocity);

        //VolumeDistibutionGPU(ref delta_heights_final, in delta_heights, in final_delta_sink, in final_delta_buld, obj_info.velocity);


        //ErosionAlgorithm(ref delta_heights_final, in delta_heights, in nodes_coordinates);

        for (int i = 1; i < terrain_size - 1; i++)
        {
            for (int j = 1; j < terrain_size - 1; j++)
            {
                heights_ter[j, i] += delta_heights_final[j, i];
            }
        }

        ErosionAlgorithmGPU(ref heights_ter, in delta_heights, num_erosion_iter);

        coef_correctness = 0.0f;
        for (int i = 1; i < terrain_size - 1; i++)
        {
            for (int j = 1; j < terrain_size - 1; j++)
            {
                coef_correctness += (heights_ter[j, i] - default_height[j, i]);
            }
        }

        //needUpdate = true;
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

        for (int i = 0; i < terrain_size; i++)
        {
            for (int j = 0; j < terrain_size; j++)
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
        for (int i = 1; i < terrain_size - 1; i++)
        {
            for (int j = 1; j < terrain_size - 1; j++)
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
        for (int i = 1; i < terrain_size - 1; i++)
        {
            for (int j = 1; j < terrain_size - 1; j++)
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
        for (int i = 1; i < terrain_size - 1; i++)
        {
            for (int j = 1; j < terrain_size - 1; j++)
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
        float [,] pressure = new float[terrain_size, terrain_size];
        for (int i = 0; i < terrain_size; i++)
        {
            for (int j = 0; j < terrain_size; j++)
            {
                pressure[j, i] = (Kc * (contur / (2 * ploshad)) + Kf) * Mathf.Pow(delta_heights[j, i] * terrain.terrainData.heightmapScale.y, n);
            }
        }
        return pressure;
    }

    float[,] calc_centrality_distr(float[,] delta_heights, int n_contacts)
    {
        float summ_distance = 0.0f;
        for (int i = 0; i < terrain_size; i++)
        {
            for (int j = 0; j < terrain_size; j++)
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

        float[,] centrality_coef = new float[terrain_size, terrain_size];

        for (int i = 0; i < terrain_size; i++)
        {
            for (int j = 0; j < terrain_size; j++)
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
        float[,] final_pressure = new float[terrain_size, terrain_size];
        for (int i = 0; i < terrain_size; i++)
        {
            for (int j = 0; j < terrain_size; j++)
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

    float[,] DetectIntersection(out List<Vector3> nodes_coordinates, float[,] heights_ter, ObjectInfo objinfo)
    {
        nodes_coordinates = new List<Vector3>();

        float[,] delta_heights = new float[terrain_size, terrain_size];
        Vector3 start = terrain.GetComponent<Transform>().position;
        Vector3 end = terrain.GetComponent<Transform>().position + new Vector3(terrain.terrainData.size.x, 0, terrain.terrainData.size.z);

        float step_x = (float)((end.x - start.x) / (terrain_size - 1));
        float step_z = (float)((end.z - start.z) / (terrain_size - 1));

        int from_1, from_2;
        int to_1, to_2;

        from_1 = (int)((objinfo.position.x - objinfo.size.x / 2)/step_x);
        from_2 = (int)((objinfo.position.z - objinfo.size.z / 2) / step_z);

        to_1 = (int)((objinfo.position.x + objinfo.size.x / 2) / step_x);
        to_2 = (int)((objinfo.position.z + objinfo.size.z / 2) / step_z);

        for (int i = from_1; i < to_1; i++)
        {
            for (int j = from_2; j < to_2; j++)
            {
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

    void BaseDeform(ref float[,] delta_heights_final, out float[,] final_delta_sink, out float[,] final_delta_buld, Vector3 Velocity,in float[,] delta_heights)
    {
        final_delta_sink = new float[terrain_size, terrain_size];
        final_delta_buld = new float[terrain_size, terrain_size];

        for (int i = 0; i < terrain_size - 1; i++)
        {
            for (int j = 0; j < terrain_size - 1; j++)
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

    void TerrainAirCorrection(ref float[,] final_delta_sink, in float[,] heights_ter)
    {
        for (int i = 0; i < terrain_size - 1; i++)
        {
            for (int j = 0; j < terrain_size - 1; j++)
            {
                if (final_delta_sink[j, i] == 0)
                    continue;

                float global_delta = (default_height[j, i] - (heights_ter[j, i] - final_delta_sink[j, i])) * scale_y;

                if (global_delta > Ma || global_delta < 0)
                {
                    continue;
                }

                final_delta_sink[j, i] = 0;
            }
        }
    }

    void VolumeDistibution(ref float[,] delta_heights_final, in float[,] delta_heights, in float[,] final_delta_sink, in float[,] final_delta_buld, Vector3 Velocity)
    {
        Vector3 start = ter_position;
        Vector3 end = ter_position + new Vector3(size_x, 0, size_z);

        float step_x = (float)((end.x - start.x) / (terrain_size - 1));
        float step_z = (float)((end.z - start.z) / (terrain_size - 1));

        node_classification[,] classification = new node_classification[terrain_size, terrain_size];
        List<Vector3> borders_coordinates = new List<Vector3>();
        classification = detect_border(classification, delta_heights, ref borders_coordinates);

        int borders_count = borders_coordinates.Count();

        //Распределеяем вытесненный объем по границам, чтобы потом с помощью эрозии сгладить 
        float[] borders_volume_distr_sink = new float[borders_count];
        float[] borders_volume_distr_buld = new float[borders_count];
        double[] delta_height_border = new double[borders_count];

        for (int i = 0; i < terrain_size - 1; i++)
        {
            for (int j = 0; j < terrain_size - 1; j++)
            {
                if (delta_heights[j, i] == 0)
                    continue;

                double coef_check = 0;
                double coef_check_1 = 0;

                double[] dist = new double[borders_count];
                double[] distr_sink_buld = new double[borders_count];
                double summ_dist = 0;
                double summ_distr_sink_buld = 0;
                for (int k = 0; k < borders_count; k++)
                {
                    dist[k] = (step_z * j - step_x * borders_coordinates[k].x) * (step_z * j - step_x * borders_coordinates[k].x) +
                              (step_x * i - step_x * borders_coordinates[k].z) * (step_x * i - step_x * borders_coordinates[k].z) +
                              (delta_heights[j, i] * scale_y - borders_coordinates[k].y * scale_y) *
                              (delta_heights[j, i] * scale_y - borders_coordinates[k].y * scale_y);
                    summ_dist += dist[k];

                    Vector3 vec_dist = new Vector3(step_z * j - step_z * borders_coordinates[k].x,
                    delta_heights[j, i] * scale_y - borders_coordinates[k].y * scale_y,
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

                for (int k = 0; k < borders_count; k++)
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

        for (int k = 0; k < borders_count; k++)
        {
            delta_heights_final[(int)borders_coordinates[k].x, (int)borders_coordinates[k].z] = (float)delta_height_border[k];
        }
    }

    void ErosionAlgorithm(ref float[,] delta_heights_final, in float[,] delta_heights, in List<Vector3> node_coordinates)
    {
        //Здесь я буду делать алгоритм эрозии
        float res = (size_x / terrain_size);
        float dzlim = (res * Mathf.Tan(Mathf.Deg2Rad * fi)) / scale_y;

        float[,] delta_erosion = new float[terrain_size, terrain_size];
        for (int k = 0; k < num_erosion_iter; k++)
        {
            for (int i = matrix_erosion_size / 2; i < terrain_size - matrix_erosion_size / 2; i++)
            {
                for (int j = matrix_erosion_size / 2; j < terrain_size - matrix_erosion_size / 2; j++)
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

    void ErosionAlgorithmGPU(ref float[,] heights_ter, in float[,] delta_heights, int num_erosion_iter)
    {
        terrain_size = terrain.terrainData.heightmapResolution;

        float res = (terrain.terrainData.size.x / terrain_size);
        float dzlim = (res * Mathf.Tan(Mathf.Deg2Rad * fi)) / terrain.terrainData.heightmapScale.y;

        output = new ComputeBuffer(terrain_size * terrain_size, sizeof(float));

        //delta = new ComputeBuffer(terrain_size * terrain_size, sizeof(float));
        //delta.SetData(delta_heights);

        for (int i = 0; i < num_erosion_iter; i++)
        {
            output.SetData(heights_ter);

            int KerID = ErrosionShader.FindKernel("Errosion");
            ErrosionShader.SetBuffer(KerID, "output_terrain", output);
            //ErrosionShader.SetBuffer(KerID, "cur_delta", delta);
            ErrosionShader.SetInt("terrain_size", terrain_size);
            ErrosionShader.SetFloat("dzlim", dzlim);

            //ErrosionShader.SetFloat("size_x", objinfo.ObjectSize.x);
            //ErrosionShader.SetFloat("size_z", objinfo.ObjectSize.z);
            //ErrosionShader.SetFloat("position_x", objinfo.ObjectPosition.x);
            //ErrosionShader.SetFloat("position_z", objinfo.ObjectPosition.z);

            ErrosionShader.Dispatch(KerID, terrain_size / 8, terrain_size / 8, 1);
            output.GetData(heights_ter);
        }
    }

    struct double3
    {
        public double x;
        public double y;
        public double z;

        public double3(double _x, double _y, double _z)
        {
            x = _x;
            y = _y;
            z = _z;
        }

    }

    void VolumeDistibutionGPU(ref float[,] delta_heights_final, in float[,] delta_heights, in float[,] final_delta_sink, in float[,] final_delta_buld, Vector3 Velocity)
    {
        Vector3 start = ter_position;
        Vector3 end = ter_position + new Vector3(size_x, 0, size_z);

        float step_x = (float)((end.x - start.x) / (terrain_size - 1));
        float step_z = (float)((end.z - start.z) / (terrain_size - 1));

        node_classification[,] classification = new node_classification[terrain_size, terrain_size];
        List<double3> borders_coordinates = new List<double3>();
        detect_border_2(classification, delta_heights, ref borders_coordinates);

        int borders_count = borders_coordinates.Count();

        if (borders_count == 0)
            return;

        //Распределеяем вытесненный объем по границам, чтобы потом с помощью эрозии сгладить 
        double[] delta_height_border = new double[borders_count];
        double[] dist = new double[borders_count];
        double[] distr_sink_buld = new double[borders_count];

        ComputeBuffer output_delta_borderb = new ComputeBuffer(borders_count, sizeof(double));
        ComputeBuffer distb = new ComputeBuffer(borders_count, sizeof(double));
        ComputeBuffer distr_buldb = new ComputeBuffer(borders_count, sizeof(double));
        ComputeBuffer border_coordsb = new ComputeBuffer(borders_count, sizeof(double) *3);

        ComputeBuffer final_delta_sinkb = new ComputeBuffer(terrain_size* terrain_size, sizeof(float));
        ComputeBuffer final_delta_buldb = new ComputeBuffer(terrain_size* terrain_size, sizeof(float));

        delta = new ComputeBuffer(terrain_size * terrain_size, sizeof(float));
        delta.SetData(delta_heights);

        output_delta_borderb.SetData(delta_height_border);
        distb.SetData(dist);
        distr_buldb.SetData(distr_sink_buld);
        border_coordsb.SetData(borders_coordinates.ToArray());

        final_delta_sinkb.SetData(final_delta_sink);
        final_delta_buldb.SetData(final_delta_buld);

        int KerID = ErrosionShader.FindKernel("VolDistr");
        ErrosionShader.SetBuffer(KerID, "border_coords", border_coordsb);
        ErrosionShader.SetBuffer(KerID, "final_delta_sink", final_delta_sinkb);
        ErrosionShader.SetBuffer(KerID, "final_delta_buld", final_delta_buldb);
        ErrosionShader.SetBuffer(KerID, "dist", distb);
        ErrosionShader.SetBuffer(KerID, "distr_buld", distr_buldb);
        ErrosionShader.SetBuffer(KerID, "output_delta_border", output_delta_borderb);
        ErrosionShader.SetBuffer(KerID, "cur_delta", delta);

        ErrosionShader.SetInt("border_count", borders_count);
        ErrosionShader.SetFloat("step_x", step_x);
        ErrosionShader.SetFloat("step_z", step_z);
        ErrosionShader.SetFloat("scale_y", scale_y);
        ErrosionShader.SetFloat("velocity_x", Velocity.x);
        ErrosionShader.SetFloat("velocity_y", Velocity.y);
        ErrosionShader.SetFloat("velocity_z", Velocity.z);

        ErrosionShader.Dispatch(KerID, terrain_size / 8, terrain_size / 8, 1);
        output_delta_borderb.GetData(delta_height_border);

        for (int k = 0; k < borders_count; k++)
        {
            delta_heights_final[(int)borders_coordinates[k].x, (int)borders_coordinates[k].z] = (float)delta_height_border[k];
        }
    }

    node_classification[,] detect_border_2(node_classification[,] nodes, float[,] delta_heights, ref List<double3> borders_coordinates)
    {
        //Определили границы
        for (int i = 1; i < terrain_size - 1; i++)
        {
            for (int j = 1; j < terrain_size - 1; j++)
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
                    borders_coordinates.Add(new double3(j, delta_heights[j, i], i));
                    //heights_ter[j, i] +=  0.005f; <---- оталдочная проверка корректности найденных границ
                }
            }
        }

        return nodes;
    }

}
