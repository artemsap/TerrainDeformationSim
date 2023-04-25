using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ObjectInformation : MonoBehaviour
{
    // Start is called before the first frame update

    public Vector3 ObjectVelocity;
    public Vector3 ObjectPosition;
    public Vector3 ObjectSize;

    Vector3 prevPosition;

    void Start()
    {
        ObjectSize = GetComponent<Collider>().bounds.size;
    }

    // Update is called once per frame
    void Update()
    {
        Vector3 currentPosition = transform.position;
        float deltaTime = Time.deltaTime;

        ObjectVelocity.x = (currentPosition.x - prevPosition.x) / deltaTime;
        ObjectVelocity.y = (currentPosition.y - prevPosition.y) / deltaTime;
        ObjectVelocity.z = (currentPosition.z - prevPosition.z) / deltaTime;

        prevPosition = currentPosition;

        ObjectPosition = currentPosition;
    }
}
