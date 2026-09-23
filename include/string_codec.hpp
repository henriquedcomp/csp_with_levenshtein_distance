#pragma once

#include <string_view>
#include <string>
#include <vector>
#include <unordered_map>

using namespace std;

class StringCodec
{
private:
    inline static string_view alphabet = "";
    inline static unordered_map<char, int> encode_map;
    inline static unordered_map<int, char> decode_map;

public:
    static void set_alphabet(string_view a)
    {
        alphabet = a;

        for (int i = 0; i < a.length(); i++)
        {
            encode_map[a[i]] = i;
            decode_map[i] = a[i];
        }
    }

    void encode(string &str, vector<int> &result)
    {
        for (int i = 0; i < str.size(); i++)
            result.push_back(encode_map[str[i]]);
    }

    void decode(vector<int> &str, string &result)
    {
        for (int i = 0; i < str.size(); i++)
            result.push_back(decode_map[str[i]]);
    }
};
