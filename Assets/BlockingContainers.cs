using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace Assembly_CSharp
{
    public class BlockingQueue<T>
    {
        Queue<T> que = new Queue<T>();
        Semaphore sem = new Semaphore(0, Int32.MaxValue);

        public void Enqueue(T item)
        {
            lock (que)
            {
                que.Enqueue(item);
            }

            sem.Release();
        }

        public T Dequeue()
        {
            sem.WaitOne();

            lock (que)
            {
                return que.Dequeue();
            }
        }

        public int Count()
        {
            return que.Count();
        }
        public T Peek()
        {
            return que.Peek();
        }
    }
}
