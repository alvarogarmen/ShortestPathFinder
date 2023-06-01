#ifndef UNTITLED_APQ_SAVING_HH
#define UNTITLED_APQ__SAVINGHH

#include <vector>
#include <unordered_map>
#include <iostream>

class APQSaving {
public:
    void insertNode(double node, double priority);
    void shiftDown(int i);
    void shiftUp(int i);
    std::pair<double, double> getMin();
    double popMin();
    void decreaseKey(double node, double value);
    bool contains(double node);
    double getPrev(double node);  // New function to retrieve the previous node

    int size() { return heap.size(); }
    bool isEmpty() { return heap.empty(); }

    std::vector<std::pair<double, double>> heap;
    std::unordered_map<double, int> index;
    std::unordered_map<double, double> prev;  // New map to store the previous node
};

void APQSaving::insertNode(double node, double priority) {
    heap.emplace_back(node, priority);
    index[node] = heap.size() - 1;
    shiftUp(heap.size() - 1);
    prev[node] = -1;  // Initialize the prev value to -1
}

void APQSaving::shiftDown(int i) {
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

void APQSaving::shiftUp(int i) {
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

double APQSaving::popMin() {
    if (isEmpty()) {
        throw std::out_of_range("Priority queue is empty.");
    }
    double minValue = heap[0].first;
    index.erase(minValue);
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

std::pair<double, double> APQSaving::getMin() {
    if (isEmpty()) {
        throw std::out_of_range("Priority queue is empty.");
    }
    return heap[0];
}

void APQSaving::decreaseKey(double node, double value) {
    if (index.find(node) == index.end()) {
        throw std::invalid_argument("Value not in priority queue.");
    }
    int i = index[node];
    if (value > heap[i].second) {
        throw std::invalid_argument("New priority must be less than or equal to current priority.");
    }
    heap[i].second = value;
    shiftUp(i);

    // Update the prev value of the node
    double parent = (i - 1) / 2;
    prev[node] = (parent >= 0) ? heap[parent].first : -1;
}

bool APQSaving::contains(double node) {
    return index.find(node) != index.end();
}

double APQSaving::getPrev(double node) {
    if (prev.find(node) == prev.end()) {
        throw std::invalid_argument("Node not found in the priority queue.");
    }
    return prev[node];
}

#endif
