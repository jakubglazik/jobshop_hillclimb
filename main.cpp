#include <iostream>
#include <vector>
#include <utility>
#include <stack>
#include <fstream>
#include <array>
#include <map>
#include <string>
#include <time.h>
#include <algorithm>
#include <unordered_set>
#include <chrono>
#include <thread>

using namespace std;
using namespace chrono;

int jobs, machines, all_jobs;

class WeightedGraph {
public:
    int vertices;
    vector<vector<array<int,3>>> adjacencyList;  // pair<destination, weight>
    vector<int> longestPath;

    WeightedGraph(int V) {
        vertices = V;
        adjacencyList.resize(V);
        for(int i=0;i<V;i++){
            longestPath.push_back(-1);
        }
    }

    void addEdge(int src, int dest, int weight,int type) {
        adjacencyList[src].push_back({dest, weight,type});
    }

    void topologicalSortUtil(int v, vector<bool>& visited, stack<int>& stack) {
        //DFS - algorytm Tarjana
        visited[v] = true;

        for (const auto& neighbor : adjacencyList[v]) {
            int nextVertex = neighbor[0];
            if (!visited[nextVertex]) {
                topologicalSortUtil(nextVertex, visited, stack);
            }
        }
        
        stack.push(v);
    }

    int topologicalSort() {
        vector<bool> visited(vertices, false);
        stack<int> s;

        for (int i = 0; i < vertices; ++i) { //do usuniecia
            if (!visited[i]) {
                topologicalSortUtil(i, visited, s);
            }
        }

        longestPath[s.top()] = 0;

        while (!s.empty()) {
            int u = s.top();
            //cout << u << endl;
            s.pop();

            if (longestPath[u] != -1) {
                for (const auto& neighbor : adjacencyList[u]) {
                    int v = neighbor[0];
                    int weight = neighbor[1];
                    if (longestPath[v] < longestPath[u] + weight) {
                            longestPath[v] = longestPath[u] + weight;
                    }
                }
            }
        }

        return longestPath[jobs*machines+1];
    
    }
        bool isCyclicUtil(int v, vector<bool>& visited, unordered_set<int>& recursionStack) {
            visited[v] = true;
            recursionStack.insert(v);

            for (const auto& adj_neighbor : adjacencyList[v]) {
                int neighbor = adj_neighbor[0];
                if (!visited[neighbor]) {
                    if (isCyclicUtil(neighbor, visited, recursionStack)) {
                        return true;
                    }
                } else if (recursionStack.find(neighbor) != recursionStack.end()) { //cykl
                    return true;
                }
            }

            recursionStack.erase(v);
            return false;
        }

    bool isCyclic() {
        vector<bool> visited(vertices, false);
        unordered_set<int> recursionStack;
        for (int i = 0; i < vertices; ++i) {
            if (!visited[i]) {
                if (isCyclicUtil(i, visited, recursionStack)) {
                    return true;
                }
            }
        }
        return false;
    }

};

    void SwapVertexes(map<int, vector<int>> &machine_map){
    vector<int> tmpDisjunVertex;
    int randomMachine = rand() % (machines);
    tmpDisjunVertex = machine_map[randomMachine];

    int index1 = rand() % tmpDisjunVertex.size();
    int randomTask1 = tmpDisjunVertex[index1];
    tmpDisjunVertex.erase(tmpDisjunVertex.begin() + index1);

    int index2 = rand() % tmpDisjunVertex.size();
    int randomTask2 = tmpDisjunVertex[index2];
    tmpDisjunVertex.erase(tmpDisjunVertex.begin() + index2);
    
    tmpDisjunVertex.clear();
    //cout<<"Random tasks: "<<randomTask1<<" "<<randomTask2<<endl;

    for(int i=0;i<machine_map[randomMachine].size();i++){
        if(machine_map[randomMachine][i] == randomTask1){
            machine_map[randomMachine][i] = randomTask2;
        }else if(machine_map[randomMachine][i] == randomTask2){
            machine_map[randomMachine][i] = randomTask1;
        }
    }
}

void PrintList(WeightedGraph &Graph){
    for(int vertex=0;vertex<machines*jobs+2;vertex++){
        cout<<vertex<<"."<<endl;
        for(int edge=0;edge<Graph.adjacencyList[vertex].size();edge++){
            cout<<Graph.adjacencyList[vertex][edge][0]<<" "<<Graph.adjacencyList[vertex][edge][1]<<" "<<Graph.adjacencyList[vertex][edge][2]<<endl;
        }
    }

}

void initDisjunctiveEgdes(WeightedGraph &Graph, map<int, vector<int>> &machine_map, vector<int> durations){
    for(int i=0;i<machines;i++){
        for(int j=0;j<jobs-1;j++){
            Graph.addEdge(machine_map[i][j],machine_map[i][j+1],durations[machine_map[i][j]],1);

        }
    }
}

void read_tailard(ifstream& inputFile, WeightedGraph &Graph, map<int, vector<int>> &machine_map, vector<int> &durations){
    int temp_time, temp_machine;

    string tmp;
    
    getline(inputFile, tmp);
    getline(inputFile, tmp);

    durations.push_back(0);

    for (int i = 0; i < jobs; ++i) {
        int shiftID = i*(machines);
        Graph.addEdge(0, shiftID+1, 0,0);
        for (int j=shiftID+1;j<machines+shiftID;j++) {
            inputFile >> temp_time;
            durations.push_back(temp_time);
            Graph.addEdge(j, j+1, temp_time,0);
        }
        inputFile >> temp_time;
        Graph.addEdge(machines + shiftID, (jobs)*(machines) + 1, temp_time,0);
        durations.push_back(temp_time);
    }
    durations.push_back(0);

    if(all_jobs!=jobs){
        getline(inputFile, tmp);
        for(int i = 0; i < all_jobs-jobs; i++){
            getline(inputFile, tmp);
        }
    }

    inputFile>>tmp;
    
    for (int i = 0; i < jobs; ++i) {
        int shiftID = i*(machines);
        for (int j=shiftID+1;j<machines+shiftID;j++) {
            inputFile >> temp_machine;
            machine_map[temp_machine-1].push_back(j);
        }
        inputFile >> temp_machine;
        machine_map[temp_machine-1].push_back(machines+shiftID);
    }

    

    inputFile.close();
}

void read_orlib(ifstream& inputFile, WeightedGraph &Graph, map<int, vector<int>> &machine_map, vector<int> &durations){
    int machine_n, duration;

    durations.push_back(0);
    for (int i = 0; i < jobs; ++i) {
        int shiftID = i*(machines);
        Graph.addEdge(0, shiftID+1, 0,0);
        for(int j=shiftID+1;j<machines+shiftID;j++)
        {
            inputFile >> machine_n >> duration;

            machine_map[machine_n].push_back(j);
            Graph.addEdge(j, j+1, duration,0);
            durations.push_back(duration);
        }
        inputFile >> machine_n >> duration;
        machine_map[machine_n].push_back(machines+shiftID);
        Graph.addEdge(machines + shiftID, (jobs)*(machines) + 1, duration,0);
        durations.push_back(duration);
    }
    
    durations.push_back(0);
    inputFile.close();

}

void czytajMape(const map<int, vector<int>>& mapa) {
    // Iteracja przez elementy mapy
    for (const auto &entry : mapa) {
        cout << "Klucz: " << entry.first << ", Wartości: ";
        
        // Iteracja przez elementy wektora
        for (const auto &value : entry.second) {
            cout << value << " ";
        }
        
        cout << endl;
    }
}

bool compareTasks(int task1, int task2) {
    return (task1-1) % machines < (task2-1) % machines;
}

void sort_machine_map(map<int, vector<int>> &mapa){
    for(int i=0;i<machines;i++)
    {
        sort(mapa[i].begin(), mapa[i].end(), compareTasks);
    }
    
}

void toFile(WeightedGraph &Graph){
    ofstream outFile("output.txt");
    outFile<<Graph.longestPath[jobs*machines+1]<<endl;

    for(int i=1;i<=jobs*machines;i++){
        outFile<<Graph.longestPath[i]<<" ";
        if(i%machines==0) outFile<<endl;
    }
    outFile.close();
}   

int main(int argc, char** argv) {
    srand(time(NULL));
    if(argc != 4)
    {
        cout << "Zla ilosc parametrow!";
        return 1;
    }  
    int nr_jobs = atoi(argv[3]);
    int iterations = 10000;
    bool bad_neighbors;


    ifstream inputFile(argv[2]);

    if (!inputFile.is_open()) {
        cout << "Nie mozna otworzyc pliku!" << endl;
        exit(1);
    }


    inputFile >> jobs >> machines;

    all_jobs = jobs;

    if(nr_jobs>jobs){
        cout<<"Zbyt mala ilosc zadan w pliku wejsciowym";
        inputFile.close();
        return 1;
    }
    else if(nr_jobs > 0) jobs = nr_jobs;
    else if(nr_jobs < -1)
    {
        cout << "Ilosc elementow musi byc >0";
        inputFile.close();
        return 1;
    }
    
    WeightedGraph baseGraph(jobs*machines+2);
    WeightedGraph mainGraph(jobs*machines+2);
    WeightedGraph neighborGraph(jobs*machines+2);
    WeightedGraph tmpNeighborGraph(jobs*machines+2);
    
    map<int, vector<int>> tmp_machine_map;
    map<int, vector<int>> machine_map;
    map<int, vector<int>> neighbor_machine_map;

    vector<int> durations; //source and sink
    int mainSolution,tmpSolution,neighborSolution;

    if(atoi(argv[1]) == 0)
    {
        read_orlib(inputFile,baseGraph, machine_map, durations);
    }
    else if(atoi(argv[1]) == 1)
    {
        read_tailard(inputFile,baseGraph, machine_map, durations);
    }
    else
    {
        cout << "Niepoprawny parametr typu instancji! wybierz 0 lub 1" << endl;
        exit(-1);
    }


    if(jobs==1){
        baseGraph.topologicalSort();
        toFile(baseGraph);
        return 0;
    }

    mainGraph = baseGraph;

    sort_machine_map(machine_map);
    
    initDisjunctiveEgdes(mainGraph,machine_map,durations);

    mainSolution = mainGraph.topologicalSort();

    neighborGraph = mainGraph;
    neighbor_machine_map = machine_map;
    neighborSolution = mainSolution;
    
    auto start_time = high_resolution_clock::now();
    cout<<mainSolution<<endl;
    while(true){
        bad_neighbors = true;
        for(int i=0;i<iterations;i++){
            tmpNeighborGraph=baseGraph;
            tmp_machine_map = machine_map;
            SwapVertexes(tmp_machine_map);
            initDisjunctiveEgdes(tmpNeighborGraph,tmp_machine_map,durations);
            if(!tmpNeighborGraph.isCyclic()){
                tmpSolution = tmpNeighborGraph.topologicalSort();
                if(neighborSolution>tmpSolution){
                    neighborSolution = tmpSolution;
                    neighborGraph = tmpNeighborGraph;
                    neighbor_machine_map = tmp_machine_map;
                    bad_neighbors = false;
                }
            }

            auto current_time = high_resolution_clock::now();
            auto elapsed_time = duration_cast<seconds>(current_time - start_time).count();
                
            if(elapsed_time >= 180) {
                cout << "Przekroczono limit czasu: "<<elapsed_time/60<<" min"<<endl;
                toFile(neighborGraph);
                return 0;
            }
            //cout<<elapsed_time<<endl;
            //cout<<bad_neighbors<<" "<<iterations<<endl;
            // if(bad_neighbors<=iterations*0.4){
            //     cout << "Nie znaleziono lepszego sasiada"<< endl;
            //     break;
            // }
        }
        mainGraph = neighborGraph;
        machine_map = neighbor_machine_map;

        if(bad_neighbors){
            cout<<"Nie znaleziono lepszych sasiadow danego rozwiazania"<<endl;
            toFile(mainGraph);
            return 0;
        }
    }
}