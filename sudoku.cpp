#include <stdio.h>
#include <cmath>
#include <cstdint>
#include <random>
#include <array>
#include <numeric>
#include <algorithm>
#include <stack>
#include <cstring>
#include <sstream>
#include <optional>
#include <tuple>
#include <map>

constexpr int N = 9;
constexpr int Ncar = sqrt(N);
constexpr int Tcar = Ncar;
constexpr int N2 = N * N;
constexpr int Ind = 18;
int largeur = 0;

typedef uint_fast8_t case_t;

typedef case_t vect[N];

using Grid = std::array<std::array<case_t, N>, N>;

char *vide;

template <typename T>
void afficher(std::array<std::array<T, N>, N> &grille)
{
    int w = (N + Ncar - 1) * (largeur + 1) - 1;
    char buf[w + 1];
    std::fill_n(buf, w, '-');
    buf[w] = 0;

    for (int y = 0; y < N; y++)
    {
        if (y > 0 && y % Tcar == 0)
            printf("%s\n", buf);

        for (int x = 0; x < N; x++)
        {
            if (x > 0 && x % Tcar == 0)
                printf("| ");

            T val = grille[y][x];

            if (val == 0)
                printf("%s ", vide);
            else
                printf("%*d ", largeur, val);
        }

        printf("\n");
    }
}

bool contient(vect v, case_t val, int i)
{
    for (int k = 0; k < i; k++)
    {
        if (v[k] == val)
            return true;
    }

    return false;
}

bool verif_ligne(Grid &grille, int y)
{
    vect exist;
    int i = 0;

    for (int x = 0; x < N; x++)
    {
        case_t val = grille[y][x];

        if (val == 0)
            continue;

        if (contient(exist, val, i))
        {
            //printf("échec ligne (%d présent deux fois sur %d)\n", val, y);
            return true;
        }

        exist[i++] = val;
    }

    return false;
}

bool verif_col(Grid &grille, int x)
{
    vect exist;
    int i = 0;

    for (int y = 0; y < N; y++)
    {
        case_t val = grille[y][x];

        if (val == 0)
            continue;

        if (contient(exist, val, i))
        {
            //printf("échec colonne (%d présent deux fois sur %d)\n", val, y);
            return true;
        }

        exist[i++] = val;
    }

    return false;
}

bool verif_carre(Grid &grille, int x, int y)
{
    vect exist; // @TODO: si Tcar² != N ?
    int i = 0;

    const int cx = (x / Tcar) * Tcar;
    const int cy = (y / Tcar) * Tcar;

    for (int y = 0; y < Tcar; y++)
        for (int x = 0; x < Tcar; x++)
        {
            case_t val = grille[cy + y][cx + x];

            if (val == 0)
                continue;

            if (contient(exist, val, i))
            {
                //printf("échec carré (%d présent deux fois sur %d,%d)\n", val, cx, cy);
                return true;
            }

            exist[i++] = val;
        }

    return false;
}

inline bool verif(Grid &grille, int x, int y)
{
    return verif_ligne(grille, y) || verif_col(grille, x) || verif_carre(grille, x, y);
}

case_t *get_lin(Grid &grille)
{
    return &grille[0][0];
}

void generer(Grid &grille)
{
    static std::random_device rd;
    static std::mt19937 mt(rd());
    static std::uniform_int_distribution<int> dist(0, N2 - 1);

    std::fill_n(get_lin(grille), N2, 0);
    auto grille_lin = get_lin(grille);

    std::vector<case_t> vals(N);
    std::iota(vals.begin(), vals.end(), 1);

    for (int k = 0; k < Ind;)
    {
        int i;

        do
        {
            i = dist(mt);
        } while (grille_lin[i] != 0);

        int l = i / N;
        int c = i % N;
        std::shuffle(vals.begin(), vals.end(), mt);

        for (case_t n : vals)
        {
            grille_lin[i] = n;

            if (!verif(grille, c, l))
            {
                k++;
                break;
            }
        }
    }
}

static int numberOfSetBits(uint32_t i)
{
    i = i - ((i >> 1) & 0x55555555);
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
    return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

static unsigned rightone32(uint32_t n)
{
    const uint32_t a = 0x05f66a47; /* magic number, found by brute force */
    static const unsigned decode[32] = {0, 1, 2, 26, 23, 3, 15, 27, 24, 21, 19, 4, 12, 16, 28, 6, 31, 25, 22, 14, 20, 18, 11, 5, 30, 13, 17, 10, 29, 9, 8, 7};
    n = a * (n & (-n));
    return decode[n >> 27];
}

class Solver
{
public:
    Solver()
    {
        zero_fill();
    }

    Solver(Grid &grille, std::optional<std::vector<int>> killer = std::nullopt)
    {
        if (this->killer = killer.has_value())
        {
            zero_fill();
            this->grille_killer = grille;
            this->killer_vals = killer.value();
            this->killer_sums = new killer_val[this->killer_vals.size()];
            printf("calcul des cellules");
            for (int i = 0; i < killer_vals.size(); i++)
            {
                std::vector<cell_pos> cells;

                for (int y = 0; y < N; y++)
                    for (int x = 0; x < N; x++)
                    {
                        if (grille_killer[y][x] == i)
                            cells.push_back((cell_pos){x, y});
                    }

                killer_cells.push_back(cells);
            }
            printf(" OK\n");
            /*this->grille[0][0] = 2;
            this->grille[1][0] = 3;
            this->grille[2][0] = 7;
            this->grille[3][0] = 5;
            this->grille[4][0] = 1;
            this->grille[5][0] = 9;*/
            //this->grille[0] = {2, 1, 5, 6, 4, 7, 3, 9, 8};
            //this->grille[1] = {3, 6, 8, 9, 5, 2, 1, 7, 4};
            //this->grille[2] = {7,9,4,3,8,1,6,5,2};
        }
        else
        {
            this->grille = grille;
        }
    }

    ~Solver()
    {
        if (killer)
            delete[] killer_sums;
    }

    Grid &get_grille()
    {
        return this->grille;
    }

    void solve()
    {
        if (killer)
        {
            killer_pre();

            solve_elim();

            solve_exhaust(0);
        }
        else
        {
            solve_elim();

            uint32_t cset = 0;

            for (int y = 0; y < N; y++)
                for (int x = 0; x < N; x++)
                    if (grille[y][x] != 0)
                        cset++;

            if (cset < N2)
            {
                solve_exhaust(cset);
            }
        }
    }

private:
    struct killer_val
    {
        int sum;
        int count;
    };

    struct cell_pos
    {
        int x;
        int y;
    };

    Grid grille;
    Grid grille_killer;
    bool killer;
    std::vector<int> killer_vals;
    std::vector<std::vector<cell_pos>> killer_cells;
    killer_val *killer_sums;
    uint32_t mask = ((1 << N) - 1) << 1;
    std::array<std::array<uint32_t, N>, N> excl = {0};
    uint32_t *excl_lin = &excl[0][0];

    void zero_fill()
    {
        std::fill_n(&this->grille[0][0], N2, 0);
    }

    void exclude(int x, int y, int val)
    {
        int bs = 1 << val;
        for (int i = 0; i < N; i++)
        {
            excl[i][x] |= bs;
            excl[y][i] |= bs;
        }

        const int cx = (x / Tcar) * Tcar;
        const int cy = (y / Tcar) * Tcar;

        for (int ly = 0; ly < Tcar; ly++)
            for (int lx = 0; lx < Tcar; lx++)
                excl[cy + ly][cx + lx] |= bs;
    }

    void exclude_line(int y, int val)
    {
        int bs = 1 << val;
        for (int i = 0; i < N; i++)
        {
            excl[y][i] |= bs;
        }
    }

    void exclude_col(int x, int val)
    {
        int bs = 1 << val;
        for (int i = 0; i < N; i++)
        {
            excl[i][x] |= bs;
        }
    }

    void include(int x, int y, int val)
    {
        excl[y][x] &= ~(1 << val);
    }

    void fill_in(int x, int y, int val, const char *msg)
    {
        grille[y][x] = val;
        exclude(x, y, val);
        printf("%-28s -> %d en (%d, %d)\n", msg, val, x, y);
    }

    void solve_elim()
    {
        printf("résolution par élimination\n");

        for (int y = 0; y < N; y++)
            for (int x = 0; x < N; x++)
            {
                int val = grille[y][x];
                if (val == 0)
                    continue;
                exclude(x, y, val);
            }

        bool updated;
        uint32_t pre_count = 0;
        do
        {
            updated = false;
            pre_count++;

            for (int y = 0; y < N; y++)
                for (int x = 0; x < N; x++)
                {
                    if (grille[y][x] != 0)
                        continue;

                    uint32_t val = excl[y][x];

                    if (val == mask)
                    {
                        afficher(grille);
                        printf("CASE IMPOSSIBLE (%d, %d)\n", x, y);
                        exit(1);
                    }

                    val = (~val) & mask;

                    if (numberOfSetBits(val) == 1)
                    {
                        fill_in(x, y, rightone32(val), "élimination exclusion");
                        updated = true;
                    }
                }

            for (int y = 0; y < N; y++)
            {
                for (int i = 1; i <= N; i++)
                {
                    int occ = 0, pos = -1;
                    for (int x = 0; x < N; x++)
                    {
                        if (grille[y][x] != 0)
                            continue;

                        if (!(excl[y][x] & (1 << i)))
                        {
                            if (++occ > 1)
                                goto multiple_l;
                            pos = x;
                        }
                    }

                    if (occ == 1)
                    {
                        fill_in(pos, y, i, "élimination ligne");
                        updated = true;
                    }

                multiple_l:
                    continue;
                }
            }

            for (int x = 0; x < N; x++)
            {
                for (int i = 1; i <= N; i++)
                {
                    int occ = 0, pos = -1;
                    for (int y = 0; y < N; y++)
                    {
                        if (grille[y][x] != 0)
                            continue;

                        if (!(excl[y][x] & (1 << i)))
                        {
                            if (++occ > 1)
                                goto multiple_c;
                            pos = y;
                        }
                    }

                    if (occ == 1)
                    {
                        fill_in(x, pos, i, "élimination colonne");
                        updated = true;
                    }

                multiple_c:
                    continue;
                }
            }

            for (int cy = 0; cy < Ncar; cy++)
                for (int cx = 0; cx < Ncar; cx++)
                {
                    int cpy = cy * Tcar;
                    int cpx = cx * Tcar;

                    for (int i = 1; i <= N; i++)
                    {
                        int occ = 0, posx = -1, posy = -1;

                        for (int y = 0; y < 3; y++)
                            for (int x = 0; x < 3; x++)
                            {
                                int ry = cpy + y;
                                int rx = cpx + x;

                                if (grille[ry][rx] != 0)
                                    continue;

                                if (!(excl[ry][rx] & (1 << i)))
                                {
                                    if (++occ > 1)
                                        goto multiple_s;

                                    posy = ry;
                                    posx = rx;
                                }
                            }

                        if (occ == 1)
                        {
                            fill_in(posx, posy, i, "élimination carré");
                            updated = true;
                        }

                    multiple_s:
                        struct
                        {
                            int x;
                            int y;
                        } occs[9];
                        occ = -1;

                        bool all_line = true;
                        bool all_col = true;

                        for (int y = 0; y < 3; y++)
                            for (int x = 0; x < 3; x++)
                            {
                                int ry = cpy + y;
                                int rx = cpy + x;

                                if (grille[ry][rx] != 0)
                                    continue;

                                if (!(excl[ry][rx] & (1 << i)))
                                {
                                    if (occ >= 0)
                                    {
                                        // ligne
                                        if (all_line && ry != occs[occ].y)
                                            all_line = false;

                                        // colonne
                                        if (all_col && rx != occs[occ].x)
                                            all_col = false;

                                        if (!all_line && !all_col)
                                            goto none_same;
                                    }

                                    occs[++occ] = {rx, ry};
                                }
                            }

                    three_occ:
                        if (occ > 0 && occ <= 2)
                        {
                            if (all_line)
                            {
                                printf("exclusion alignement ligne   de %d en %d et ", i, occs[0].y);
                                exclude_line(occs[0].y, i);
                                printf("réinclusion en ");
                                for (int j = 0; j <= occ; j++)
                                {
                                    printf("(%d, %d)  ", occs[j].x, occs[j].y);
                                    include(occs[j].x, occs[j].y, i);
                                }
                                printf("\n");
                            }
                            else if (all_col)
                            {
                                printf("exclusion alignement colonne de %d en %d et ", i, occs[0].x);
                                exclude_col(occs[0].x, i);
                                printf("réinclusion en ");
                                for (int j = 0; j <= occ; j++)
                                {
                                    printf("(%d, %d)  ", occs[j].x, occs[j].y);
                                    include(occs[j].x, occs[j].y, i);
                                }
                                printf("\n");
                            }
                        }

                    none_same:
                        continue;
                    }
                }

        } while (updated);

        printf("nb pre:   %10d\n", pre_count);
    }

    bool check_killer_regions()
    {
        std::fill_n(killer_sums, killer_vals.size(), (killer_val){0, 0});

        for (int y = 0; y < N; y++)
            for (int x = 0; x < N; x++)
            {
                int pos = grille_killer[y][x];
                int val = grille[y][x];
                int target = killer_vals[pos];

                if (val != 0)
                {
                    killer_sums[pos].sum += val;
                    killer_sums[pos].count++;
                }

                if (killer_sums[pos].sum > target)
                    return true;

                if (killer_sums[pos].count == killer_cells[pos].size() && killer_sums[pos].sum != target)
                    return true;
            }

        return false;
    }

    void killer_pre()
    {
        printf("prétraitement killer\n");
        typedef std::map<int, std::vector<int>> pmap;
        std::vector<std::map<int, std::vector<std::vector<int>>>> possibles =
            {
                {{3, {{1, 2}}},
                 {4, {{1, 3}}},
                 {5, {{1, 4}, {2, 3}}},
                 {6, {{1, 5}, {2, 4}}},
                 {7, {{1, 6}, {2, 5}, {3, 4}}},
                 {8, {{1, 7}, {2, 6}, {3, 5}}},
                 {9, {{1, 8}, {2, 7}, {3, 6}, {4, 5}}},
                 {10, {{1, 9}, {2, 8}, {3, 7}, {4, 6}}},
                 {11, {{2, 9}, {3, 8}, {4, 7}, {5, 6}}},
                 {12, {{3, 9}, {4, 8}, {5, 7}}},
                 {13, {{4, 9}, {5, 8}, {6, 7}}},
                 {14, {
                          {5, 9},
                          {6, 8},
                      }},
                 {15, {{6, 9}, {7, 8}}},
                 {16, {{7, 9}}},
                 {17, {{8, 9}}}},
                {{6, {{1, 2, 3}}}, {7, {{1, 2, 4}}}, {8, {{1, 2, 5}, {1, 3, 4}}}, {9, {{1, 2, 6}, {1, 3, 5}, {2, 3, 4}}}, {10, {{1, 2, 7}, {1, 3, 6}, {1, 4, 5}, {2, 3, 5}}}, {11, {{1, 2, 8}, {1, 3, 7}, {1, 4, 6}, {2, 3, 6}, {2, 4, 5}}}, {12, {{1, 2, 9}, {1, 3, 8}, {1, 4, 7}, {1, 5, 6}, {2, 3, 7}, {2, 4, 6}, {3, 4, 5}}}, {13, {{1, 3, 9}, {1, 4, 8}, {1, 5, 7}, {2, 3, 8}, {2, 4, 7}, {2, 5, 6}, {3, 4, 6}}}, {14, {{1, 4, 9}, {1, 5, 8}, {1, 6, 7}, {2, 3, 9}, {2, 4, 8}, {2, 5, 7}, {3, 4, 7}, {3, 5, 6}}}, {15, {{1, 5, 9}, {1, 6, 8}, {2, 4, 9}, {2, 5, 8}, {2, 6, 7}, {3, 4, 8}, {3, 5, 7}, {4, 5, 6}}}, {16, {{1, 6, 9}, {1, 7, 8}, {2, 5, 9}, {2, 6, 8}, {3, 4, 9}, {3, 5, 8}, {3, 6, 7}, {4, 5, 7}}}, {17, {{1, 7, 9}, {2, 6, 9}, {2, 7, 8}, {3, 5, 9}, {3, 6, 8}, {4, 5, 8}, {4, 6, 7}}}, {18, {{1, 8, 9}, {2, 7, 9}, {3, 6, 9}, {3, 7, 8}, {4, 5, 9}, {4, 6, 8}, {5, 6, 7}}}, {19, {{2, 8, 9}, {3, 7, 9}, {4, 6, 9}, {4, 7, 8}, {5, 6, 8}}}, {20, {{3, 8, 9}, {4, 7, 9}, {5, 6, 9}, {5, 7, 8}}}, {21, {{4, 8, 9}, {5, 7, 9}, {6, 7, 8}}}, {22, {{5, 8, 9}, {6, 7, 9}}}, {23, {{6, 8, 9}}}, {24, {{7, 8, 9}}}},
                {{10, {{1, 2, 3, 4}}}, {11, {{1, 2, 3, 5}}}, {12, {{1, 2, 3, 6}, {1, 2, 4, 5}}}, {13, {{1, 2, 3, 7}, {1, 2, 4, 6}, {1, 3, 4, 5}}}, {14, {{1, 2, 3, 8}, {1, 2, 4, 7}, {1, 2, 5, 6}, {1, 3, 4, 6}, {2, 3, 4, 5}}}, {15, {{1, 2, 3, 9}, {1, 2, 4, 8}, {1, 2, 5, 7}, {1, 3, 4, 7}, {1, 3, 5, 6}, {2, 3, 4, 6}}}, {16, {{1, 2, 4, 9}, {1, 2, 5, 8}, {1, 2, 6, 7}, {1, 3, 4, 8}, {1, 3, 5, 7}, {1, 4, 5, 6}, {2, 3, 4, 7}, {2, 3, 5, 6}}}, {17, {{1, 2, 5, 9}, {1, 2, 6, 8}, {1, 3, 4, 9}, {1, 3, 5, 8}, {1, 3, 6, 7}, {1, 4, 5, 7}, {2, 3, 4, 8}, {2, 3, 5, 7}, {2, 4, 5, 6}}}, {18, {{1, 2, 6, 9}, {1, 2, 7, 8}, {1, 3, 5, 9}, {1, 3, 6, 8}, {1, 4, 5, 8}, {1, 4, 6, 7}, {2, 3, 4, 9}, {2, 3, 5, 8}, {2, 3, 6, 7}, {2, 4, 5, 7}, {3, 4, 5, 6}}}, {19, {{1, 2, 7, 9}, {1, 3, 6, 9}, {1, 3, 7, 8}, {1, 4, 5, 9}, {1, 4, 6, 8}, {1, 5, 6, 7}, {2, 3, 5, 9}, {2, 3, 6, 8}, {2, 4, 5, 8}, {2, 4, 6, 7}, {3, 4, 5, 7}}}, {20, {{1, 2, 8, 9}, {1, 3, 7, 9}, {1, 4, 6, 9}, {1, 4, 7, 8}, {1, 5, 6, 8}, {2, 3, 6, 9}, {2, 3, 7, 8}, {2, 4, 5, 9}, {2, 4, 6, 8}, {2, 5, 6, 7}, {3, 4, 5, 8}, {3, 4, 6, 7}}}, {21, {{1, 3, 8, 9}, {1, 4, 7, 9}, {1, 5, 6, 9}, {1, 5, 7, 8}, {2, 3, 7, 9}, {2, 4, 6, 9}, {2, 4, 7, 8}, {2, 5, 6, 8}, {3, 4, 5, 9}, {3, 4, 6, 8}, {3, 5, 6, 7}}}, {22, {{1, 4, 8, 9}, {1, 5, 7, 9}, {1, 6, 7, 8}, {2, 3, 8, 9}, {2, 4, 7, 9}, {2, 5, 6, 9}, {2, 5, 7, 8}, {3, 4, 6, 9}, {3, 4, 7, 8}, {3, 5, 6, 8}, {4, 5, 6, 7}}}, {23, {{1, 5, 8, 9}, {1, 6, 7, 9}, {2, 4, 8, 9}, {2, 5, 7, 9}, {2, 6, 7, 8}, {3, 4, 7, 9}, {3, 5, 6, 9}, {3, 5, 7, 8}, {4, 5, 6, 8}}}, {24, {{1, 6, 8, 9}, {2, 5, 8, 9}, {2, 6, 7, 9}, {3, 4, 8, 9}, {3, 5, 7, 9}, {3, 6, 7, 8}, {4, 5, 6, 9}, {4, 5, 7, 8}}}, {25, {{1, 7, 8, 9}, {2, 6, 8, 9}, {3, 5, 8, 9}, {3, 6, 7, 9}, {4, 5, 7, 9}, {4, 6, 7, 8}}}, {26, {{2, 7, 8, 9}, {3, 6, 8, 9}, {4, 5, 8, 9}, {4, 6, 7, 9}, {5, 6, 7, 8}}}, {27, {{3, 7, 8, 9}, {4, 6, 8, 9}, {5, 6, 7, 9}}}, {28, {{4, 7, 8, 9}, {5, 6, 8, 9}}}, {29, {{5, 7, 8, 9}}}, {30, {{6, 7, 8, 9}}}},
                {
                    {15, {{1, 2, 3, 4, 5}}},
                    {16, {{1, 2, 3, 4, 6}}},
                    {17, {{1, 2, 3, 4, 7}, {1, 2, 3, 5, 6}}},
                    {18, {{1, 2, 3, 4, 8}, {1, 2, 3, 5, 7}, {1, 2, 4, 5, 6}}},
                    {19, {{1, 2, 3, 4, 9}, {1, 2, 3, 5, 8}, {1, 2, 3, 6, 7}, {1, 2, 4, 5, 7}, {1, 3, 4, 5, 6}}},
                    {20, {{1, 2, 3, 5, 9}, {1, 2, 3, 6, 8}, {1, 2, 4, 5, 8}, {1, 2, 4, 6, 7}, {1, 3, 4, 5, 7}, {2, 3, 4, 5, 6}}},
                    {21, {{1, 2, 3, 6, 9}, {1, 2, 3, 7, 8}, {1, 2, 4, 5, 9}, {1, 2, 4, 6, 8}, {1, 2, 5, 6, 7}, {1, 3, 4, 5, 8}, {1, 3, 4, 6, 7}, {2, 3, 4, 5, 7}}},
                    {22, {{1, 2, 3, 7, 9}, {1, 2, 4, 6, 9}, {1, 2, 4, 7, 8}, {1, 2, 5, 6, 8}, {1, 3, 4, 5, 9}, {1, 3, 4, 6, 8}, {1, 3, 5, 6, 7}, {2, 3, 4, 5, 8}, {2, 3, 4, 6, 7}}},
                    {23, {{1, 2, 3, 8, 9}, {1, 2, 4, 7, 9}, {1, 2, 5, 6, 9}, {1, 2, 5, 7, 8}, {1, 3, 4, 6, 9}, {1, 3, 4, 7, 8}, {1, 3, 5, 6, 8}, {1, 4, 5, 6, 7}, {2, 3, 4, 5, 9}, {2, 3, 4, 6, 8}, {2, 3, 5, 6, 7}}},
                    {24, {{1, 2, 4, 8, 9}, {1, 2, 5, 7, 9}, {1, 2, 6, 7, 8}, {1, 3, 4, 7, 9}, {1, 3, 5, 6, 9}, {1, 3, 5, 7, 8}, {1, 4, 5, 6, 8}, {2, 3, 4, 6, 9}, {2, 3, 4, 7, 8}, {2, 3, 5, 6, 8}, {2, 4, 5, 6, 7}}},
                    {25, {{1, 2, 5, 8, 9}, {1, 2, 6, 7, 9}, {1, 3, 4, 8, 9}, {1, 3, 5, 7, 9}, {1, 3, 6, 7, 8}, {1, 4, 5, 6, 9}, {1, 4, 5, 7, 8}, {2, 3, 4, 7, 9}, {2, 3, 5, 6, 9}, {2, 3, 5, 7, 8}, {2, 4, 5, 6, 8}, {3, 4, 5, 6, 7}}},
                    {26, {{1, 2, 6, 8, 9}, {1, 3, 5, 8, 9}, {1, 3, 6, 7, 9}, {1, 4, 5, 7, 9}, {1, 4, 6, 7, 8}, {2, 3, 4, 8, 9}, {2, 3, 5, 7, 9}, {2, 3, 6, 7, 8}, {2, 4, 5, 6, 9}, {2, 4, 5, 7, 8}, {3, 4, 5, 6, 8}}},
                    {27, {{1, 2, 7, 8, 9}, {1, 3, 6, 8, 9}, {1, 4, 5, 8, 9}, {1, 4, 6, 7, 9}, {1, 5, 6, 7, 8}, {2, 3, 5, 8, 9}, {2, 3, 6, 7, 9}, {2, 4, 5, 7, 9}, {2, 4, 6, 7, 8}, {3, 4, 5, 6, 9}, {3, 4, 5, 7, 8}}},
                    {28, {{1, 3, 7, 8, 9}, {1, 4, 6, 8, 9}, {1, 5, 6, 7, 9}, {2, 3, 6, 8, 9}, {2, 4, 5, 8, 9}, {2, 4, 6, 7, 9}, {2, 5, 6, 7, 8}, {3, 4, 5, 7, 9}, {3, 4, 6, 7, 8}}},
                    {29, {{1, 4, 7, 8, 9}, {1, 5, 6, 8, 9}, {2, 3, 7, 8, 9}, {2, 4, 6, 8, 9}, {2, 5, 6, 7, 9}, {3, 4, 5, 8, 9}, {3, 4, 6, 7, 9}, {3, 5, 6, 7, 8}}},
                    {30, {{1, 5, 7, 8, 9}, {2, 4, 7, 8, 9}, {2, 5, 6, 8, 9}, {3, 4, 6, 8, 9}, {3, 5, 6, 7, 9}, {4, 5, 6, 7, 8}}},
                    {31, {{1, 6, 7, 8, 9}, {2, 5, 7, 8, 9}, {3, 4, 7, 8, 9}, {3, 5, 6, 8, 9}, {4, 5, 6, 7, 9}}},
                    {32, {{2, 6, 7, 8, 9}, {3, 5, 7, 8, 9}, {4, 5, 6, 8, 9}}},
                    {33, {{3, 6, 7, 8, 9}, {4, 5, 7, 8, 9}}},
                    {34, {{4, 6, 7, 8, 9}}},
                    {35, {{5, 6, 7, 8, 9}}},
                },
                {
                    {21, {{1, 2, 3, 4, 5, 6}}},
                    {22, {{1, 2, 3, 4, 5, 7}}},
                    {23, {{1, 2, 3, 4, 5, 8}, {1, 2, 3, 4, 6, 7}}},
                    {24, {{1, 2, 3, 4, 5, 9}, {1, 2, 3, 4, 6, 8}, {1, 2, 3, 5, 6, 7}}},
                    {25, {{1, 2, 3, 4, 6, 9}, {1, 2, 3, 4, 7, 8}, {1, 2, 3, 5, 6, 8}, {1, 2, 4, 5, 6, 7}}},
                    {26, {{1, 2, 3, 4, 7, 9}, {1, 2, 3, 5, 6, 9}, {1, 2, 3, 5, 7, 8}, {1, 2, 4, 5, 6, 8}, {1, 3, 4, 5, 6, 7}}},
                    {27, {{1, 2, 3, 4, 8, 9}, {1, 2, 3, 5, 7, 9}, {1, 2, 3, 6, 7, 8}, {1, 2, 4, 5, 6, 9}, {1, 2, 4, 5, 7, 8}, {1, 3, 4, 5, 6, 8}, {2, 3, 4, 5, 6, 7}}},
                    {28, {{1, 2, 3, 5, 8, 9}, {1, 2, 3, 6, 7, 9}, {1, 2, 4, 5, 7, 9}, {1, 2, 4, 6, 7, 8}, {1, 3, 4, 5, 6, 9}, {1, 3, 4, 5, 7, 8}, {2, 3, 4, 5, 6, 8}}},
                    {29, {{1, 2, 3, 6, 8, 9}, {1, 2, 4, 5, 8, 9}, {1, 2, 4, 6, 7, 9}, {1, 2, 5, 6, 7, 8}, {1, 3, 4, 5, 7, 9}, {1, 3, 4, 6, 7, 8}, {2, 3, 4, 5, 6, 9}, {2, 3, 4, 5, 7, 8}}},
                    {30, {{1, 2, 3, 7, 8, 9}, {1, 2, 4, 6, 8, 9}, {1, 2, 5, 6, 7, 9}, {1, 3, 4, 5, 8, 9}, {1, 3, 4, 6, 7, 9}, {1, 3, 5, 6, 7, 8}, {2, 3, 4, 5, 7, 9}, {2, 3, 4, 6, 7, 8}}},
                    {31, {{1, 2, 4, 7, 8, 9}, {1, 2, 5, 6, 8, 9}, {1, 3, 4, 6, 8, 9}, {1, 3, 5, 6, 7, 9}, {1, 4, 5, 6, 7, 8}, {2, 3, 4, 5, 8, 9}, {2, 3, 4, 6, 7, 9}, {2, 3, 5, 6, 7, 8}}},
                    {32, {{1, 2, 5, 7, 8, 9}, {1, 3, 4, 7, 8, 9}, {1, 3, 5, 6, 8, 9}, {1, 4, 5, 6, 7, 9}, {2, 3, 4, 6, 8, 9}, {2, 3, 5, 6, 7, 9}, {2, 4, 5, 6, 7, 8}}},
                    {33, {{1, 2, 6, 7, 8, 9}, {1, 3, 5, 7, 8, 9}, {1, 4, 5, 6, 8, 9}, {2, 3, 4, 7, 8, 9}, {2, 3, 5, 6, 8, 9}, {2, 4, 5, 6, 7, 9}, {3, 4, 5, 6, 7, 8}}},
                    {34, {{1, 3, 6, 7, 8, 9}, {1, 4, 5, 7, 8, 9}, {2, 3, 5, 7, 8, 9}, {2, 4, 5, 6, 8, 9}, {3, 4, 5, 6, 7, 9}}},
                    {35, {{1, 4, 6, 7, 8, 9}, {2, 3, 6, 7, 8, 9}, {2, 4, 5, 7, 8, 9}, {3, 4, 5, 6, 8, 9}}},
                    {36, {{1, 5, 6, 7, 8, 9}, {2, 4, 6, 7, 8, 9}, {3, 4, 5, 7, 8, 9}}},
                    {37, {{2, 5, 6, 7, 8, 9}, {3, 4, 6, 7, 8, 9}}},
                    {38, {{3, 5, 6, 7, 8, 9}}},
                    {39, {{4, 5, 6, 7, 8, 9}}},
                },
                {
                    {28, {{1, 2, 3, 4, 5, 6, 7}}},
                    {29, {{1, 2, 3, 4, 5, 6, 8}}},
                    {30, {{1, 2, 3, 4, 5, 6, 9}, {1, 2, 3, 4, 5, 7, 8}}},
                    {31, {{1, 2, 3, 4, 5, 7, 9}, {1, 2, 3, 4, 6, 7, 8}}},
                    {32, {{1, 2, 3, 4, 5, 8, 9}, {1, 2, 3, 4, 6, 7, 9}, {1, 2, 3, 5, 6, 7, 8}}},
                    {33, {{1, 2, 3, 4, 6, 8, 9}, {1, 2, 3, 5, 6, 7, 9}, {1, 2, 4, 5, 6, 7, 8}}},
                    {34, {{1, 2, 3, 4, 7, 8, 9}, {1, 2, 3, 5, 6, 8, 9}, {1, 2, 4, 5, 6, 7, 9}, {1, 3, 4, 5, 6, 7, 8}}},
                    {35, {{1, 2, 3, 5, 7, 8, 9}, {1, 2, 4, 5, 6, 8, 9}, {1, 3, 4, 5, 6, 7, 9}, {2, 3, 4, 5, 6, 7, 8}}},
                    {36, {{1, 2, 3, 6, 7, 8, 9}, {1, 2, 4, 5, 7, 8, 9}, {1, 3, 4, 5, 6, 8, 9}, {2, 3, 4, 5, 6, 7, 9}}},
                    {37, {{1, 2, 4, 6, 7, 8, 9}, {1, 3, 4, 5, 7, 8, 9}, {2, 3, 4, 5, 6, 8, 9}}},
                    {38, {{1, 2, 5, 6, 7, 8, 9}, {1, 3, 4, 6, 7, 8, 9}, {2, 3, 4, 5, 7, 8, 9}}},
                    {39, {{1, 3, 5, 6, 7, 8, 9}, {2, 3, 4, 6, 7, 8, 9}}},
                    {40, {{1, 4, 5, 6, 7, 8, 9}, {2, 3, 5, 6, 7, 8, 9}}},
                    {41, {{2, 4, 5, 6, 7, 8, 9}}},
                    {42, {{3, 4, 5, 6, 7, 8, 9}}},
                },
                {
                    {36, {{1, 2, 3, 4, 5, 6, 7, 8}}},
                    {37, {{1, 2, 3, 4, 5, 6, 7, 9}}},
                    {38, {{1, 2, 3, 4, 5, 6, 8, 9}}},
                    {39, {{1, 2, 3, 4, 5, 7, 8, 9}}},
                    {40, {{1, 2, 3, 4, 6, 7, 8, 9}}},
                    {41, {{1, 2, 3, 5, 6, 7, 8, 9}}},
                    {42, {{1, 2, 4, 5, 6, 7, 8, 9}}},
                    {43, {{1, 3, 4, 5, 6, 7, 8, 9}}},
                    {44, {{2, 3, 4, 5, 6, 7, 8, 9}}},
                },
                {{45, {{1, 2, 3, 4, 5, 6, 7, 8, 9}}}}};

        for (int i = 0; i < killer_vals.size(); i++)
        {
            if (killer_cells[i].size() == 1)
            {
                printf("région %2d triviale avec %d\n", i, killer_vals[i]);
                auto cell = killer_cells[i][0];
                grille[cell.y][cell.x] = killer_vals[i];
                continue;
            }
            int sizeIdx = killer_cells[i].size() - 2;

            if (sizeIdx < 0 || sizeIdx >= possibles.size())
                continue;

            auto map = possibles[sizeIdx];

            auto iter = map.find(killer_vals[i]);
            if (iter == map.end())
                continue;
            std::vector<std::vector<int>> vals = iter->second;
            printf("région n°%2d = %2d : ", i, killer_vals[i]);
            int mask = ~0;
            for (std::vector<int> vals2 : vals)
            {
                printf("(");
                for (int val : vals2)
                {
                    printf("%d ", val);
                    mask &= ~(1 << val);
                }
                printf("\b) ");
            }
            printf("en ");
            for (auto cell : killer_cells[i])
            {
                printf("(%d, %d) ", cell.x, cell.y);
                excl[cell.y][cell.x] = mask;
            }
            printf("\n");
        }
    }

    void solve_exhaust(int cset)
    {
        if (!killer)
            afficher(grille);

        printf("résolution exhaustive itérative\n");

        auto grille_lin = get_lin(grille);

        uint16_t *log = new uint16_t[10'000'000];
        uint32_t log_pos = 0;

        constexpr int lprog = 100;
        int ni = N2 - cset;
        char buf[lprog + 1];
        buf[lprog] = 0;

        std::stack<int> i;

        i.push(0);

        uint32_t nb_hypot = 0, nb_excl = 0;
        bool bp = false;
        while (true)
        {
            int k;

            log[log_pos++] = i.size();

            int maxs = (i.size() - 1) * lprog / ni;

            for (k = 0; k < maxs; k++)
                buf[k] = '@';
            for (; k < lprog; k++)
                buf[k] = '.';

            i.push(0);

            printf("\r%s", buf);
            fflush(stdout);

            while (i.top() < N2 && grille_lin[i.top()] != 0)
                i.top()++;

            std::stack<int> icp = i;

            if (i.top() == N2)
                break;

            while (true)
            {
                if (i.empty())
                {
                    printf("\n");
                    afficher(grille);
                    printf("PAS DE SOLUTION ! [%d]\n", nb_hypot);
                    exit(1);
                }

                case_t &act = grille_lin[i.top()];

                if (act == N)
                {
                    // rebouclage
                    act = 0;
                    i.pop();
                    continue;
                }

                act++;

                nb_hypot++;

                /*if (grille[0][0] == 2 && grille[0][2] == 5 && grille[0][3] == 6 && grille[0][4] == 4)
                    bp = true;*/

                if (bp)
                {
                    printf("\n");
                    afficher(grille);
                    getchar();
                }

                if (excl_lin[i.top()] & (1 << act))
                {
                    if (bp)
                        printf("exclu %d en (%d, %d)\n", excl_lin[i.top()], i.top() % N, i.top() / N);
                    nb_excl++;
                    continue;
                }

                int l = i.top() / N;
                int c = i.top() % N;

                if (verif(grille, c, l) || (killer && check_killer_regions()))
                {
                    // backtrack
                    continue;
                }
                else
                {
                    // réussite
                    break;
                }
            }
        }

        printf("\n");
        printf("nb hypot: %10d\n", nb_hypot);
        printf("nb excl:  %10d (%3d%%)\n", nb_excl, nb_excl * 100 / nb_hypot);

        auto file = fopen("prog.csv", "w");
        auto lptr = &log[0];
        for (uint32_t i = 0; i < log_pos; i++, lptr++)
            fprintf(file, "%d\n", *lptr);
        fclose(file);

        delete[] log;
    }
};

int get_num(char c)
{
    if (c >= '0' && c <= '9')
        return c - '0';
    if (c >= 'a' && c <= 'z')
        return 10 + c - 'a';

    printf("wtf %c\n", c);
    exit(1);
}

int main(int argc, char *argv[])
{
    static_assert(Ncar * Ncar == N, "ERREUR: Ncar incorrect");

    int n = N;
    do
    {
        largeur++;
        n /= 10;
    } while (n);

    vide = new char[largeur + 1];

    for (int i = 0; i < largeur; i++)
        vide[i] = '.';

    vide[largeur] = 0;

    Grid grille;

    auto args = &argv[1];
    bool killer = false;
    std::vector<int> killer_vals;

    if (argc > 1)
    {
        if (strcmp(args[0], "-k") == 0)
        {
            killer = true;
            args++;
            argc--;
        }

        if (argc >= 10)
        {
            std::fill_n(get_lin(grille), N2, 0);

            char **ligne = args;
            for (int i = 1; i < 10; i++, ligne++)
                for (int j = 0; j < 9; j++)
                {
                    if (killer)
                    {
                        grille[i - 1][j] = get_num((*ligne)[j]);
                    }
                    else
                    {
                        if ((*ligne)[j] >= '0' && (*ligne)[j] <= '9')
                            grille[i - 1][j] = (*ligne)[j] - '0';
                    }
                }

            if (killer)
            {
                std::stringstream in(*ligne);
                int temp;
                char ch;
                while (in >> temp)
                {
                    printf("%d ", temp);
                    killer_vals.push_back(temp);
                    in >> ch;
                }
                printf("\n");
            }
        }
    }
    else
    {
        generer(grille);
    }

    if (!killer)
        afficher(grille);

    Solver s(grille, killer ? std::optional(killer_vals) : std::nullopt);
    s.solve();

    printf("\n");

    afficher(s.get_grille());
}