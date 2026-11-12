#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <random>
#include <format>

using namespace std;
namespace fs = filesystem;

int main()
{
    // controls if an instance generated now will replace content from a file with the same name
    const bool REPLACE_INSTANCES = false;

    // controls if instances/random folder will be completely deleted and instances recreated from zero
    const bool DELETE_AND_RECREATE_INSTANCES = true;

    const int INSTANCES_PER_COMBINATION = 3;

    const string INSTANCE_PATH = "./instances";
    const string RANDOM_INSTANCE_PATH = INSTANCE_PATH + "/random";
    const string RANDOM_BINARY_INSTANCE_PATH = RANDOM_INSTANCE_PATH + "/binary";
    const string RANDOM_DNA_INSTANCE_PATH = RANDOM_INSTANCE_PATH + "/dna";
    const string RANDOM_PROTEIN_INSTANCE_PATH = RANDOM_INSTANCE_PATH + "/protein";

    const vector<int> UB_STRING_SIZE = {10, 15, 20};
    const vector<int> NUMBER_OF_STRINGS = {10, 20, 30, 40};

    const string_view BINARY_ALPHABET = "01";
    const string_view NUCLEOTIDE_ALPHABET = "ACGT";
    const string_view AMINOACID_ALPHABET = "ARNDCQEGHILKMFPSTWYV";

    if (fs::exists(RANDOM_INSTANCE_PATH) && DELETE_AND_RECREATE_INSTANCES)
        fs::remove_all(RANDOM_INSTANCE_PATH);

    random_device rd;
    mt19937 gen(rd());

    cout << "Starting BINARY random instance generation..." << endl;

    uniform_int_distribution<int>
        b_alphabet_char_dist(0, BINARY_ALPHABET.size() - 1);
    for (int x : UB_STRING_SIZE)
    {
        uniform_int_distribution<int> str_len_dist(floor(x / 2.0), x);

        for (int y : NUMBER_OF_STRINGS)
        {
            for (int i = 0; i < INSTANCES_PER_COMBINATION; i++)
            {
                const string PATH = RANDOM_BINARY_INSTANCE_PATH + format("/BINARY_{}_{}_{}.in", x, y, i);

                if (fs::exists(PATH) && !REPLACE_INSTANCES)
                {
                    cout << "Ignoring instance replace: " << PATH << endl;

                    continue;
                }

                fs::path fs_path(PATH);
                fs::path directories = fs_path.parent_path();

                if (!directories.empty())
                    fs::create_directories(directories);

                ofstream file(PATH);

                file << BINARY_ALPHABET << endl;
                file << y << endl;

                for (int j = 0; j < y; j++)
                {
                    int str_size = str_len_dist(gen);

                    string str;

                    for (int k = 0; k < str_size; k++)
                        str.push_back(BINARY_ALPHABET[b_alphabet_char_dist(gen)]);

                    file << str << endl;
                }

                file.close();
            }
        }
    }

    cout << "Starting DNA random instance generation..." << endl;

    uniform_int_distribution<int>
        n_alphabet_char_dist(0, NUCLEOTIDE_ALPHABET.size() - 1);
    for (int x : UB_STRING_SIZE)
    {
        uniform_int_distribution<int> str_len_dist(floor(x / 2.0), x);

        for (int y : NUMBER_OF_STRINGS)
        {
            for (int i = 0; i < INSTANCES_PER_COMBINATION; i++)
            {
                const string PATH = RANDOM_DNA_INSTANCE_PATH + format("/DNA_{}_{}_{}.in", x, y, i);

                if (fs::exists(PATH) && !REPLACE_INSTANCES)
                {
                    cout << "Ignoring instance replace: " << PATH << endl;

                    continue;
                }

                fs::path fs_path(PATH);
                fs::path directories = fs_path.parent_path();

                if (!directories.empty())
                    fs::create_directories(directories);

                ofstream file(PATH);

                file << NUCLEOTIDE_ALPHABET << endl;
                file << y << endl;

                for (int j = 0; j < y; j++)
                {
                    int str_size = str_len_dist(gen);

                    string str;

                    for (int k = 0; k < str_size; k++)
                        str.push_back(NUCLEOTIDE_ALPHABET[n_alphabet_char_dist(gen)]);

                    file << str << endl;
                }

                file.close();
            }
        }
    }

    cout << "Starting PROTEIN random instance generation..." << endl;

    uniform_int_distribution<int> a_alphabet_char_dist(0, AMINOACID_ALPHABET.size() - 1);
    for (int x : UB_STRING_SIZE)
    {
        uniform_int_distribution<int> str_len_dist(floor(x / 2.0), x);

        for (int y : NUMBER_OF_STRINGS)
        {
            for (int i = 0; i < INSTANCES_PER_COMBINATION; i++)
            {
                const string PATH = RANDOM_PROTEIN_INSTANCE_PATH + format("/PROTEIN_{}_{}_{}.in", x, y, i);

                if (fs::exists(PATH) && !REPLACE_INSTANCES)
                {
                    cout << "Ignoring instance replace: " << PATH << endl;

                    continue;
                }

                fs::path fs_path(PATH);
                fs::path directories = fs_path.parent_path();

                if (!directories.empty())
                    fs::create_directories(directories);

                ofstream file(PATH);

                file << AMINOACID_ALPHABET << endl;
                file << y << endl;

                for (int j = 0; j < y; j++)
                {
                    int str_size = str_len_dist(gen);

                    string str;

                    for (int k = 0; k < str_size; k++)
                        str.push_back(AMINOACID_ALPHABET[a_alphabet_char_dist(gen)]);

                    file << str << endl;
                }

                file.close();
            }
        }
    }

    cout << "Instance generation finished successfully!" << endl;

    return 0;
}
