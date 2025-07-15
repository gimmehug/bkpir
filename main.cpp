// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

#include "main.h"

using namespace std;
using namespace seal;

size_t g_N = 8192;
size_t g_num_n = 128;
size_t g_num_m = 8192;
size_t g_max_keyword_bytelength = 10;
size_t g_max_item_bytelength = 2;
size_t g_k = 4;

int main(int argc, char* argv[])
{
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            show_help(argv[0]);
            return 0;
        } else if (arg == "-N" && i + 1 < argc) {
            g_N = std::stoul(argv[++i]);
        } else if (arg == "-n" && i + 1 < argc) {
            g_num_n = std::stoul(argv[++i]);
        } else if (arg == "-m" && i + 1 < argc) {
            g_num_m = std::stoul(argv[++i]);
        } else if (arg == "-kl" && i + 1 < argc) {
            g_max_keyword_bytelength = std::stoul(argv[++i]);
        } else if (arg == "-il" && i + 1 < argc) {
            g_max_item_bytelength = std::stoul(argv[++i]);
        } else if (arg == "-k" && i + 1 < argc) {
            g_k = std::stoul(argv[++i]);
        } else {
            std::cerr << "Unknown option or missing argument: " << arg << std::endl;
            show_help(argv[0]);
            return 1;
        }
    }

    while (true)
    {
        cout << "+-----------------------------------------------+" << endl;
        cout << "| Examples             | Source Files           |" << endl;
        cout << "+----------------------+------------------------+" << endl;
        cout << "| 1. MMKPIRL           | mmkpirl.cpp            |" << endl;
        cout << "| 2. BKPIRL-Single     | bkpirl_single.cpp      |" << endl;
        cout << "| 3. BKPIRL-AND        | bkpirl_and.cpp         |" << endl;
        cout << "| 4. BKPIRL-OR         | bkpirl_or.cpp          |" << endl;
        cout << "| 5. BKPIRL-NOT        | bkpirl_not.cpp         |" << endl;
        cout << "+----------------------+------------------------+" << endl;

        int selection = 0;
        bool valid = true;
        do
        {
            cout << endl << "> Run example (1 ~ 5) or exit (0): ";
            if (!(cin >> selection))
            {
                valid = false;
            }
            else if (selection < 0 || selection > 5)
            {
                valid = false;
            }
            else
            {
                valid = true;
            }
            if (!valid)
            {
                cout << "  [Beep~~] valid option: type 0 ~ 5" << endl;
                cin.clear();
                cin.ignore(numeric_limits<streamsize>::max(), '\n');
            }
        } while (!valid);

        switch (selection)
        {
        case 1:
            mmkpirl(g_N, g_num_n, g_num_m, g_max_keyword_bytelength, 
                          g_max_item_bytelength, g_k);
            break;

        case 2:
            bkpirl_single(g_N, g_num_n, g_num_m, g_max_keyword_bytelength, 
                          g_max_item_bytelength, g_k);
            break;

        case 3:
            bkpirl_and(g_N, g_num_n, g_num_m, g_max_keyword_bytelength, 
                          g_max_item_bytelength, g_k);
            break;

        case 4:
            bkpirl_or(g_N, g_num_n, g_num_m, g_max_keyword_bytelength, 
                          g_max_item_bytelength, g_k);
            break;

        case 5:
            bkpirl_not(g_N, g_num_n, g_num_m, g_max_keyword_bytelength, 
                          g_max_item_bytelength, g_k);
            break;

        case 0:
            return 0;
        }
    }

    return 0;
}

void show_help(const char* program_name) {
    std::cout << "Usage: " << program_name << " [options]\n"
              << "Options:\n"
              << "  -h, --help         Show this help message\n"
              << "  -N  <value>        Polynomial modulus (default: 8192)\n"
              << "  -n  <value>        Number of keywords (default: 128)\n"
              << "  -m  <value>        Number of items (default: 8192)\n"
              << "  -kl <value>        Max keyword byte length (default: 10)\n"
              << "  -il <value>        Max item byte length (default: 2)\n"
              << "  -k  <value>        CW code Hamming weight (default: 4)\n";
}