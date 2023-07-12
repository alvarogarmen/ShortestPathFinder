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
    APQ(double size);
    void insertNode(double node, double priority);      //Insert a node
    void shiftDown(int i);                              //Shift down operation
    void shiftUp(int i);                                //Shift up operation
    std::pair<double, double> getMin();                 //Return the value of the node with the smallest priority
    double popMin();                                    //DeleteMin operation
    void decreaseKey(double node, double value);        //Normal decreaseKey operation
    bool contains(double node);                         //Check if element already in it

    int size(){return heap.size();}
    bool isEmpty(){return heap.empty();}


    // Vector with the node and their priority, aka our priority queue. Left is nodeId and right is the priority
    std::vector<std::pair<double, double>> heap;
    // A hash table that maps values to their indices in the heap. Only important here
    std::vector<int> index;
};

APQ::APQ(double size) {
    this->index.resize(size, -1);
}
void APQ::insertNode(double node, double priority){
    heap.emplace_back(node, priority);
    index[node] = heap.size()-1;
    shiftUp(index[node]);
}

void APQ::shiftDown(int i) {
    int size = heap.size();
    while (2 * i + 1 < size) {
        int left = 2 * i + 1;
        int right = 2 * i + 2;
        int smallest = left;
        if (right < size && heap[right].second < heap[left].second) {
            smallest = right;
        }
        if (heap[i].second > heap[smallest].second) {
            std::swap(heap[i], heap[smallest]);
            index[heap[i].first] = i;
            index[heap[smallest].first] = smallest;
            i = smallest;
        } else {
            break;
        }
    }
}

void APQ::shiftUp(int i) {
    while (i > 0) {
        int parent = (i - 1) / 2;
        if (heap[i].second < heap[parent].second) {
            std::swap(heap[i], heap[parent]);
            index[heap[i].first] = i;
            index[heap[parent].first] = parent;
            i = parent;
        } else {
            break;
        }
    }
}

double APQ::popMin() {
    if (isEmpty()) {
        throw std::out_of_range("Priority queue is empty.");
    }
    double minValue = heap[0].first;
    index[minValue] = -1;
    if (heap.size() > 1) {
        heap[0] = std::move(heap.back());
        index[heap[0].first] = 0;
        heap.pop_back();
        shiftDown(0);
    } else {
        heap.pop_back();
    }

    return minValue;
}

std::pair<double, double> APQ::getMin() {
    if (isEmpty()) {
        throw std::out_of_range("Priority queue is empty.");
    }
    return heap[0];
}

void APQ::decreaseKey(double node, double value) {
    if (index[node] == -1){
        throw std::invalid_argument("Value not in priority queue");
    }
    int i = index[node];
    if (value > heap[i].second) {
        throw std::invalid_argument("New priority must be less than or equal to current priority.");
    }
    heap[i].second = value;
    shiftUp(i);
}

bool APQ::contains(double node) {
    return index[node] != -1;
}
#endif
