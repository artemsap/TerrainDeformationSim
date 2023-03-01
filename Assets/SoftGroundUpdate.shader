Shader "Soft Ground Height Update"
{
    Properties
    {
        _brushSize("Brush size", Range(0, 1)) = 0.1
        _brushIntensity("Brush intensity", Range(-1, 1)) = 0
    }

        SubShader
    {
        Lighting Off
        Blend One Zero

        Pass
        {
            CGPROGRAM
            #include "UnityCustomRenderTexture.cginc"
            #pragma vertex CustomRenderTextureVertexShader
            #pragma fragment frag
            #pragma target 3.0

            float4 _DrawPosition;
            float _brushSize;
            float _brushIntensity;

            float4 frag(v2f_customrendertexture IN) : COLOR
            {
                float4 prevColor = tex2D(_SelfTexture2D, IN.localTexcoord.xy);
                float4 drawColor = smoothstep(_brushIntensity, _brushSize, distance(IN.localTexcoord.xy, _DrawPosition));
               
                float r = clamp(prevColor.r - abs(drawColor.r - 1), 0, 1);

                return float4(r, 0, 0, 0);
            }
            ENDCG
        }
    }
}