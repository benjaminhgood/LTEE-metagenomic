#ifndef OBJECT_POOL
#define OBJECT_POOL

#include<vector>
#include<algorithm>
#include<iostream>

template <class T>
class SharedObjectPool;

template <class T>
class SharedObjectPool{
    public:

        class ObjectRecord{
            public:
                T object;
                int count;
        };

        int capacity;
        std::vector<ObjectRecord> object_records;
        std::vector<int> free_list;

        SharedObjectPool(int initial_capacity = 1, T const & prototype_object = T()): capacity(initial_capacity+1), object_records(capacity, ObjectRecord{prototype_object,0}), free_list(capacity,0) {
            for(int i=0; i<capacity; ++i){
                free_list[i] = i;
            }
        };

        void reserve(int new_capacity){
            if(new_capacity > capacity){
                object_records.reserve(new_capacity);
                object_records.insert(object_records.end(), new_capacity-capacity, ObjectRecord{object_records.back().object, 0} );
                free_list.reserve(new_capacity);
                for(int i=capacity; i<new_capacity; ++i){
                    free_list.push_back(i);
                }
                capacity = new_capacity;
            }
        }

        void free(int object_id){ free_list.push_back(object_id); };
 
        class ObjectPtr{
            public:
                SharedObjectPool * pool;
                int object_id;
            
                void swap(ObjectPtr & object_ptr){ 
                    std::swap(pool, object_ptr.pool);
                    std::swap(object_id, object_ptr.object_id);
                };
                void swap(ObjectPtr && object_ptr){
                    std::swap(pool, object_ptr.pool);
                    std::swap(object_id, object_ptr.object_id);
                };


                ObjectPtr(): pool(nullptr), object_id(-1) { std::cout << "Shouldn't be doing this a lot!" << std::endl; };
                ObjectPtr(SharedObjectPool * pool, int object_id): pool(pool), object_id(object_id) { ++pool->object_records[object_id].count; };
                ObjectPtr(ObjectPtr const & o): pool(o.pool), object_id(o.object_id) { ++pool->object_records[object_id].count;};
                ObjectPtr(ObjectPtr && o): pool(nullptr), object_id(-1) {
                    swap(o);
                };

                ObjectPtr & operator=(ObjectPtr const & o){
                    swap(ObjectPtr(o));
                    return *this;
                }
                ObjectPtr & operator=(ObjectPtr && o){
                    swap(o);
                    return *this;
                }
                ~ObjectPtr(){ 
                    if(object_id != -1){
                        --pool->object_records[object_id].count;
                        if(pool->object_records[object_id].count == 0){
                            pool->free(object_id);
                        }
                    }
                };
                T* operator->(){ return &pool->object_records[object_id].object; };
                const T* operator->() const { return &pool->object_records[object_id].object;};  
                T& operator*() { return pool->object_records[object_id].object; };
                const T& operator*() const { return pool->object_records[object_id].object; };
        };

        ObjectPtr allocate(){
            if(free_list.empty()) reserve(capacity*1.5+1);
            ObjectPtr o(this,free_list.back());
            free_list.pop_back();
            return o;
        };
        
};

#endif
