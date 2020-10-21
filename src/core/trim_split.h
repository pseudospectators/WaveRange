
#include <vector>
std::string trim(const std::string& string, const char* trimCharacterList);
std::vector<std::string> split(std::string str, char* del);

// http://program.station.ez-net.jp/special/handbook/cpp/string/trim.asp
// http://www.martinbroadhurst.com/how-to-trim-a-stdstring.html
std::string trim(const std::string& string, const char* trimCharacterList = " \t\v\r\n")
{
    std::string result;
    std::string::size_type left = string.find_first_not_of(trimCharacterList);
    if (left != std::string::npos)
    {
        std::string::size_type right = string.find_last_not_of(trimCharacterList);
        result = string.substr(left, right - left + 1);
    }
    return result;
}

// split_find in https://www.sejuku.net/blog/49378
std::vector<std::string> split(std::string str, char* del) {
    int first = 0;
    int last = str.find_first_of(del);

    std::vector<std::string> result;

    while (first < str.size()) {
        std::string subStr(str, first, last - first);

        result.push_back(subStr);

        first = last + 1;
        last = str.find_first_of(del, first);

        if (last == std::string::npos) {
            last = str.size();
        }
    }

    return result;
}
