#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>


// Function to trim whitespace from a string
std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t");
    return str.substr(first, last - first + 1);
}

std::map<std::string, std::string> parseConfig(const std::string& filename) {
    std::ifstream file(filename);
    std::map<std::string, std::string> config;
    std::string line;

    if (!file.is_open()) {
        std::cerr << "Could not open the file!" << std::endl;
        return config;
    }

    while (getline(file, line)) {
        // Remove comments
        size_t commentPos = line.find('#');
        if (commentPos != std::string::npos) {
            line = line.substr(0, commentPos);
        }

        // Remove whitespace from both ends
        line = trim(line);

        // Skip empty lines
        if (line.empty()) continue;

        // Split the line into key and value
        size_t delimiterPos = line.find('=');
        if (delimiterPos != std::string::npos) {
            std::string key = trim(line.substr(0, delimiterPos));
            std::string value = trim(line.substr(delimiterPos + 1));
            config[key] = value;
        }
    }

    file.close();
    return config;
}


// Function to get an integer value by key
double getDoubleValue(const std::map<std::string, std::string>& config, const std::string& key) {
    auto it = config.find(key);
    if (it != config.end()) {
        try {
            return std::stod(it->second);
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid double value for key: " << key << std::endl;
        } catch (const std::out_of_range& e) {
            std::cerr << "Double value out of range for key: " << key << std::endl;
        }
    }
    return 0.0; // Default value
}

//read DB
// Function to read key-value pairs from a CSV file and store them in a map
std::map<int, int> readCSVToMap(const std::string &filename) {
    std::ifstream file(filename);
    std::string line;
    std::map<int, int> dataMap;

    // Check if the file is open
    if (file.is_open()) {
        while (getline(file, line)) {
            std::stringstream ss(line);
            std::string keyStr, valueStr;
            
            // Read the key and value from the line
            if (getline(ss, keyStr, ',') && getline(ss, valueStr, ',')) {
                int key = std::stoi(keyStr);  // Convert key to int
		int value = std::stoi(valueStr);  // Convert value to int
                dataMap[key] = value;
            }
        }
        file.close();
    } else {
        std::cerr << "Unable to open file";
    }

    return dataMap;
}

// Function to lookup value by key
int lookupValue(const std::map<int, int> &dataMap, int key) {
    auto it = dataMap.find(key);
    if (it != dataMap.end()) {
        return it->second;
    } else {
        return -1;
    }
}

//read Aexp
// Function to read data from file and return a map with run number and (Aexp, errAexp) pairs
// Define a struct to store all data for each run
struct RunData {
    double Aexp;
    double errAexp;
    double pol_He3;
    double pol_beam;
};

// Function to read data from file and return a map with run number and RunData
std::map<int, RunData> readRunDataFromFile(const std::string& filename) {
    std::ifstream inputFile(filename); // Open the file
    std::map<int, RunData> data; // Map to store run number and RunData

    if (!inputFile) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return data; // Return empty map if file can't be opened
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        std::istringstream iss(line);
        std::string temp;
        int runNumber;
        RunData runData;

        // Extract the run number
        iss >> temp >> runNumber;

        // Skip the "| Aexp:" part and extract Aexp
        iss >> temp >> temp >> runData.Aexp;

        // Skip the "| errAexp:" part and extract errAexp
        iss >> temp >> temp >> runData.errAexp;

        // Skip the "| pol_He3:" part and extract pol_He3
        iss >> temp >> temp >> runData.pol_He3;

        // Skip the "| pol_beam:" part and extract pol_beam
        iss >> temp >> temp >> runData.pol_beam;

        // Store in the map
        data[runNumber] = runData;
    }

    inputFile.close();
    return data; // Return the map with data
}





