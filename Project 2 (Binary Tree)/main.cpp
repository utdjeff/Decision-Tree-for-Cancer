
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <memory>
#include <string>

struct Node {
    int feature_index;
    int threshold;
    std::string classification;
    std::shared_ptr<Node> left_child;
    std::shared_ptr<Node> right_child;

    Node(int index, int thresh, std::string classif = "")
        : feature_index(index), threshold(thresh), classification(classif) {}

    bool is_leaf() const {
        return !classification.empty();
    }
};

std::shared_ptr<Node> build_decision_tree() {
    // Root node based on the Uniformity of Cell Size
    auto root = std::make_shared<Node>(2, 2); // feature_index 2 corresponds to Uniformity of Cell Size <= 2

    // If Uniformity of Cell Size <= 2
    root->left_child = std::make_shared<Node>(6, 3); // feature_index 6 corresponds to Bare Nuclei
    // If Bare Nuclei <= 3
    root->left_child->left_child = std::make_shared<Node>(0, 0, "Benign");
    // If Bare Nuclei > 3
    root->left_child->right_child = std::make_shared<Node>(1, 3); // feature_index 1 corresponds to Clump Thickness
    // If Clump Thickness <= 3
    root->left_child->right_child->left_child = std::make_shared<Node>(0, 0, "Benign");
    // If Clump Thickness > 3
    root->left_child->right_child->right_child = std::make_shared<Node>(7, 2); // feature_index 7 corresponds to Bland Chromatin
    // If Bland Chromatin <= 2
    root->left_child->right_child->right_child->left_child = std::make_shared<Node>(4, 3); // feature_index 4 corresponds to Marginal Adhesion
    // If Marginal Adhesion <= 3
    root->left_child->right_child->right_child->left_child->left_child = std::make_shared<Node>(0, 0, "Malignant");
    // If Marginal Adhesion > 3
    root->left_child->right_child->right_child->left_child->right_child = std::make_shared<Node>(0, 0, "Benign");
    // If Bland Chromatin > 2
    root->left_child->right_child->right_child->right_child = std::make_shared<Node>(0, 0, "Malignant");

    // If Uniformity of Cell Size > 2
    root->right_child = std::make_shared<Node>(3, 2); // feature_index 3 corresponds to Uniformity of Cell Shape
    // If Uniformity of Cell Shape <= 2
    root->right_child->left_child = std::make_shared<Node>(1, 5); // feature_index 1 corresponds to Clump Thickness
    // If Clump Thickness <= 5
    root->right_child->left_child->left_child = std::make_shared<Node>(0, 0, "Benign");
    // If Clump Thickness > 5
    root->right_child->left_child->right_child = std::make_shared<Node>(0, 0, "Malignant");
    // If Uniformity of Cell Shape > 2
    root->right_child->right_child = std::make_shared<Node>(2, 4); // feature_index 2 corresponds to Uniformity of Cell Size
    // If Uniformity of Cell Size <= 4
    root->right_child->right_child->left_child = std::make_shared<Node>(6, 2); // feature_index 6 corresponds to Bare Nuclei
    // If Bare Nuclei <= 2
    root->right_child->right_child->left_child->left_child = std::make_shared<Node>(4, 3); // feature_index 4 corresponds to Marginal Adhesion
    // If Marginal Adhesion <= 3
    root->right_child->right_child->left_child->left_child->left_child = std::make_shared<Node>(0, 0, "Benign");
    // If Marginal Adhesion > 3
    root->right_child->right_child->left_child->left_child->right_child = std::make_shared<Node>(0, 0, "Malignant");
    // If Bare Nuclei > 2
    root->right_child->right_child->left_child->right_child = std::make_shared<Node>(1, 6); // feature_index 1 corresponds to Clump Thickness
    // If Clump Thickness <= 6
    root->right_child->right_child->left_child->right_child->left_child = std::make_shared<Node>(2, 3); // feature_index 2 corresponds to Uniformity of Cell Size
    // If Uniformity of Cell Size <= 3
    root->right_child->right_child->left_child->right_child->left_child->left_child = std::make_shared<Node>(0, 0, "Malignant");
    // If Uniformity of Cell Size > 3
    root->right_child->right_child->left_child->right_child->left_child->right_child = std::make_shared<Node>(4, 5); // feature_index 4 corresponds to Marginal Adhesion
    // If Marginal Adhesion <= 5
    root->right_child->right_child->left_child->right_child->left_child->right_child->left_child = std::make_shared<Node>(0, 0, "Benign");
    // If Marginal Adhesion > 5
    root->right_child->right_child->left_child->right_child->left_child->right_child->right_child = std::make_shared<Node>(0, 0, "Malignant");
    // If Clump Thickness > 6
    root->right_child->right_child->left_child->right_child->right_child = std::make_shared<Node>(0, 0, "Malignant");
    // If Uniformity of Cell Size > 4
    root->right_child->right_child->right_child = std::make_shared<Node>(0, 0, "Malignant");

    return root;
}

std::string classify_tree(const std::vector<std::string>& instance, const std::shared_ptr<Node>& node, int& invalid_count) {
    // Check for "?" in the specified feature indices before classifying
    for (int i : {1, 2, 4, 5, 6, 7, 8, 9}) { // The indices to be checked
        if (instance[i] == "?") {
            invalid_count++; // Increment invalid count if a "?" is found
            return "Invalid"; // Return a string indicating this is an invalid instance
        }
    }

    if (node->is_leaf()) {
        return node->classification;
    }

    int feature_value = std::stoi(instance[node->feature_index]);
    if (feature_value <= node->threshold) {
        return classify_tree(instance, node->left_child, invalid_count);
    }
    else {
        return classify_tree(instance, node->right_child, invalid_count);
    }
}

std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

std::vector<std::vector<std::string>> read_csv(const std::string& filename) {
    std::vector<std::vector<std::string>> data;
    std::ifstream file(filename);

    std::string line;
    while (std::getline(file, line)) {
        std::vector<std::string> row = split(line, ',');
        data.push_back(row);
    }

    return data;
}

void write_csv(const std::string& filename, const std::vector<std::vector<std::string>>& data) {
    std::ofstream file(filename);

    for (const auto& row : data) {
        for (size_t i = 0; i < row.size(); i++) {
            file << row[i];
            if (i < row.size() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }
}

int main() {
    const std::string input_filename = "unformatted_data_v1.0.0.csv";
    const std::string output_filename = "results.csv";
    std::vector<std::vector<std::string>> data = read_csv(input_filename);
    std::vector<std::vector<std::string>> classified_data;

    auto decision_tree = build_decision_tree();

    int benign_count = 0;
    int malignant_count = 0;
    int invalid_count = 0;

    for (auto& instance : data) {
        std::string diagnosis = classify_tree(instance, decision_tree, invalid_count);
        instance.push_back(diagnosis);
        classified_data.push_back(instance);

        if (diagnosis == "Benign") {
            benign_count++;
        }
        else if (diagnosis == "Malignant") {
            malignant_count++;
        }
        // Invalid instances are already counted in the classify_tree function
    }

    write_csv(output_filename, classified_data);

    std::cout << "Total Patients Processed: " << classified_data.size() << std::endl;
    std::cout << "Total Benign: " << benign_count << std::endl;
    std::cout << "Total Malignant: " << malignant_count << std::endl;
    std::cout << "Total Invalid Instances: " << invalid_count << std::endl;

    return 0;
}
