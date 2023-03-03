using UnityEngine;

public class SoftGroundBrush : MonoBehaviour
{
    public CustomRenderTexture SoftGroundHeightMap;
    public Material SoftGroundMat;

    private Camera cam;
    private static readonly int DrawPosition = Shader.PropertyToID("_DrawPosition");

    // Start is called before the first frame update
    void Start()
    {
        SoftGroundMat.SetVector(DrawPosition, new Vector2(-1,-1));
        SoftGroundHeightMap.Initialize();
        cam = Camera.main;
    }

    // Update is called once per frame
    void Update()
    {
        if(Input.GetMouseButton(0))
        {
            Ray ray = cam.ScreenPointToRay(Input.mousePosition);
            if (Physics.Raycast(ray, out RaycastHit hit))
            {
                Vector2 hitTextureCoord = hit.textureCoord;
                SoftGroundMat.SetVector(DrawPosition, hitTextureCoord);
                SoftGroundHeightMap.Update();
            }
        }
        else if(Input.GetMouseButtonUp(0))
        { 
            SoftGroundMat.SetVector(DrawPosition, new Vector2(-1, -1));
            SoftGroundHeightMap.Update();
        }
       //
    }
}
