//
// Created by Alvaro on 19/04/2023.
//

#ifndef UNTITLED_APQ_HH
#define UNTITLED_APQ_HH
#include <vector>
#include <unordered_map>
#include <iostream>

class APQ{
public:

    void insertNode(double node, double priority);      //Insert a node
    void shiftDown(int i);                              //Shift down operation
    void shiftUp(int i);                                //Shift up operation
    std::pair<double, double> getMin();                 //Return the value of the node with the smallest priority
    double popMin();                                    //DeleteMin operation
    void decreaseKey(double node, double value);        //Normal decreaseKey operation
    bool contains(double node);                         //Check if element already in it

    int size(){return heap.size();}
    bool isEmpty(){return heap.empty();}

private:
    // Vector with the node and their priority, aka our priority queue. Left is nodeId and right is the priority
    std::vector<std::pair<double, double>> heap;
    // A hash table that maps values to their indices in the heap. Only important here
    std::unordered_map<double, int> index;
};
void APQ::insertNode(double node, double priority){
    this->heap.emplace_back(node, priority);
    shiftUp(size()-1);
}

void APQ::shiftDown(int i) {
    while (2 * i + 1 < size()) {
        int left = 2 * i + 1;
        int right = 2 * i + 2;
        int smallest = left;
        if (right < size() && heap[right].first < heap[left].first) {
            smallest = right;
        }
        if (heap[i].first > heap[smallest].first) {
            std::swap(heap[i], heap[smallest]);
            index[heap[i].second] = i;
            index[heap[smallest].second] = smallest;
            i = smallest;
        } else {
            break;
        }
    }
}

void APQ::shiftUp(int i) {
    while (i > 0) {
        int p = (i - 1) / 2;
        if (heap[p].first > heap[i].first) {
            std::swap(heap[p], heap[i]);
            index[heap[p].second] = p;
            index[heap[i].second] = i;
            i = p;
        } else {
            break;
        }
    }
}

double APQ::popMin() {
    double min_value = heap.front().second;
    index.erase(min_value);
    if(heap.size() > 1){
        heap.front() = std::move(heap.back());
        index[heap.front().second] = 0;
        heap.pop_back();
        shiftDown(0);
    }
    else{
        heap.pop_back();
    }
    return min_value;
}

std::pair<double, double> APQ::getMin() {
    if (heap.empty()){
        std::cout<<"Heap is empty"<<std::endl;
        return {0,0};
    }
    return heap.front();
}

void APQ::decreaseKey(double node, double value) {
    if (index.find(node) == index.end()) {
        throw std::invalid_argument("Value not in priority queue.");
    }
    int i = index[node];
    if (value > heap[i].first) {
        throw std::invalid_argument("New priority must be less than or equal to current priority.");
    }
    heap[i].first = value;
    shiftUp(i);
}

bool APQ::contains(double node) {
    for (auto & i : heap){
        if(i.first==node){
            return true;
        }
    }
    return false;
}
#endif
