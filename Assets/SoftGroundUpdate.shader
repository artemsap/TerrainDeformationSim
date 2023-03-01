Shader "Soft Ground Height Update"
{
    Properties
    {
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

            float4 frag(v2f_customrendertexture IN) : COLOR
            {
                float4 prevColor = tex2D(_SelfTexture2D, IN.localTexcoord.xy);
                float4 drawColor = smoothstep(0, 0.1, distance(IN.localTexcoord.xy, _DrawPosition));
               
                return clamp(prevColor - abs(drawColor - 1), float4(0, 0, 0, 0), float4(1, 1, 1, 1));
            }
            ENDCG
        }
    }
}