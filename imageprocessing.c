#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "imageprocessing.h"
#define _255 255
// Functia flip_horizontal va oglindi o imaginie pe orizontala
int ***flip_horizontal(int ***image, int N, int M) {
    // Cream o copia a lui image si o alocam
    int ***cp_image = (int***)malloc(N * sizeof(int**));
    for (int i = 0; i < N; i++) {
        cp_image[i] = (int**)malloc(M * sizeof(int*));
        for (int j = 0; j < M; j++) {
            cp_image[i][j] = (int*)malloc(3 * sizeof(int));
        }
    }
    if (!cp_image) {
        return NULL;
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            // iteram prin imagine si punem pe elementul cp_image[i][M - j - 1][k] elementul image[i][j][k] unde k=0,1,2
            cp_image[i][M - j -1][0] = image[i][j][0];
            cp_image[i][M - j -1][1] = image[i][j][1];
            cp_image[i][M - j -1][2] = image[i][j][2];
        }
    }
    // dezalocam image
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            free(image[i][j]);
        }
        free(image[i]);
    }
    free(image);
    return cp_image;
}
// Aceasta functie va returna imaginea primita ca parametru rotita cu 90 de grade spre stanga
int ***rotate_left(int ***image, int N, int M) {
    // alocam o copie a image pe care functia o va returna la sfarsitul executiei sale
    int ***rotate_image = (int***)malloc(M*sizeof(int**));
    for (int i = 0; i < M; i++) {
        rotate_image[i] = (int**)malloc(N*sizeof(int*));
        for (int j = 0; j < N; j++) {
            rotate_image[i][j] = (int*)malloc(3*sizeof(int));
        }
    }
    if (!rotate_image) {
        return NULL;
    }
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < 3; k++) {
                // parcurgem matricea image si elementul image[j][M - i - 1][k]
                // il atribuim copiei imaginii de pe pozitia curenta
                rotate_image[i][j][k] = image[j][M - i - 1][k];
            }
        }
    }
    // dezalocam matricea image
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            free(image[i][j]);
        }
        free(image[i]);
    }
    free(image);
    return rotate_image;
}
// Aceasta functie va returna imaginea obtinuta prin crop
// care incepe de la pozitia cu coordonatele x,y si va avea marimea h,w
int ***crop(int ***image, int N, int M, int x, int y, int h, int w) {
    // initializam o matrice in care vom mentine imaginea obtinuta prin crop a image
    int ***crop_image = (int***)malloc(h * sizeof(int**));
    for (int i = 0; i < h; i++) {
        crop_image[i] = (int**)malloc(w * sizeof(int*));
        for (int j = 0; j < w; j++) {
            crop_image[i][j] = (int*)malloc(3 * sizeof(int));
        }
    }
    if (!crop_image) {
        return NULL;
    }
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            for (int k = 0; k < 3; k++) {
                // Mai intai de a face crop la imagine verificam daca elementul
                // [i+y][j+x] nu iese in afara imaginii initiale
                if (i+y >= N || j+x >= M) {
                    continue;
                }
                // Aplicam crop asupra imaginii initiale
                crop_image[i][j][k] = image[i+y][j+x][k];
            }
        }
    }
    // dezalocam image
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            free(image[i][j]);
        }
        free(image[i]);
    }
    free(image);
    return crop_image;
}
// aceasta functie va mari imaginea primita ca paramentru linii atât deasupra cât și dedesubt
// și cu cols coloane atât la stânga cât și la dreapta
// iar pixelii nou creati trebuie sa aiba culoarea data de new_R, new_g, new_B
int ***extend(int ***image, int N, int M, int rows, int cols, int new_R, int new_G, int new_B) {
    // initializam si creeam o matrice care va contine imaginea marita
    // deoarece marim imaginea cu un numar de rows sus si jos atunci vom avea marimea matricii N + 2 * rows
    // asemanator si pentru M, M -> M + 2 * cols
    int ***extend_image = (int***)calloc((N + 2 * rows), sizeof(int**));
    int new_N = N + 2 * rows, new_M = M + 2 * cols;
    for (int i = 0; i < new_N; i++) {
        extend_image[i] = (int**)calloc(new_M, sizeof(int*));
        for (int j = 0; j < new_M; j++) {
            extend_image[i][j] = (int*)calloc(3, sizeof(int));
        }
    }
    if (!extend_image) {
        return NULL;
    }
    for (int i = rows; i < N + rows; i++) {
        for (int j = cols; j < cols + M; j++) {
            for (int k = 0; k < 3; k++) {
                // marim imaginea cu un numar anumit de linii si coloane
                extend_image[i][j][k] = image[i-rows][j-cols][k];
            }
        }
    }
    // coloram pixelii care se afla la stanga imaginii
    for (int i = 0; i < new_N; i++) {
        for (int j = 0; j < cols; j++) {
            extend_image[i][j][0] = new_R;
            extend_image[i][j][1] = new_G;
            extend_image[i][j][2] = new_B;
        }
    }
    // coloram pixelii care se afla la dreapta imaginii
    for (int i = 0; i < new_N; i++) {
        for (int j = cols + M; j< new_M; j++) {
            extend_image[i][j][0] = new_R;
            extend_image[i][j][1] = new_G;
            extend_image[i][j][2] = new_B;
        }
    }
    // coloram pixelii care se afla deasupra imaginii
    for (int i = 0; i < rows; i++) {
        for (int j = cols; j < cols + M; j++) {
            extend_image[i][j][0] = new_R;
            extend_image[i][j][1] = new_G;
            extend_image[i][j][2] = new_B;
        }
    }
    // coloram pixelii care se afla dedesubtul imaginii
    for (int i = rows + N; i < new_N; i++) {
        for (int j = cols; j < cols + M; j++) {
            extend_image[i][j][0] = new_R;
            extend_image[i][j][1] = new_G;
            extend_image[i][j][2] = new_B;
        }
    }
    // dezalocam imaginea
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            free(image[i][j]);
        }
        free(image[i]);
    }
    free(image);
    return extend_image;
}
// aceasta functie trebuie sa dea copy-paste la image_src peste image_dst
int ***paste(int ***image_dst, int N_dst, int M_dst, int ***image_src, int N_src, int M_src, int x, int y) {
     for (int i = 0; i < N_src; i++) {
        // daca ne aflam intrun index i care este invalid adica in afara image_dst atunci iesim din ciclu
        if (i >= N_dst - y) {
            break;
        }
        for (int j = 0; j < M_src; j++) {
            // daca ne aflam intrun index j care este invalid adica in afara image_dst atunci iesim din ciclu
            if (j >= M_dst-x) {
                break;
            }
            for (int k = 0; k < 3; k++) {
                // dam copy-paste image_src peste image_dst incepand de la niste coordonate x,y
                image_dst[i+y][j+x][k] = image_src[i][j][k];
            }
        }
    }
    return image_dst;
}
// aceasta functie aplica un filtru asupra unei imagini
int ***apply_filter(int ***image, int N, int M, float **filter, int filter_size) {
    // initializam o matrice care va contine imaginea asupra careia va fi aplicat filtrul
    int ***filter_image = (int***)malloc(N * sizeof(int**));
    for (int i = 0; i < N; i++) {
        filter_image[i] = (int**)malloc(M * sizeof(int*));
        for (int j = 0; j < M; j++) {
            filter_image[i][j] = (int*)malloc(3*sizeof(int));
        }
    }
    if (!filter_image) {
        return NULL;
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            float new_R = 0, new_G = 0, new_B = 0;
            for (int i1 = 0; i1 < filter_size; i1++) {
                for (int j1 = 0; j1 < filter_size; j1++) {
                    int i2 = i - filter_size/2 + i1;
                    int j2 = j - filter_size/2 + j1;
                    // verificam daca ne aflam pe un pixel care se afla in interiorul imaginii
                    if (i2 >= 0 && i2 < N && j2 >= 0 && j2 < M) {
                        // calculam noile R, G, B pentru a le putea aplica asupra imaginii
                        new_R += (float)image[i2][j2][0] * filter[i1][j1];
                        new_G += (float)image[i2][j2][1] * filter[i1][j1];
                        new_B += (float)image[i2][j2][2] * filter[i1][j1];
                    }
                }
            }
            // aplicam filtrul asupra imagnii
            filter_image[i][j][0] = (int)new_R;
            filter_image[i][j][1] = (int)new_G;
            filter_image[i][j][2] = (int)new_B;
            for (int k = 0; k < 3; k++) {
                // verificam daca filtrul aplicat asupra imaginii eeste valid deoarece culorarea variaza intre 0 si 255
                if (filter_image[i][j][k] < 0) {
                    filter_image[i][j][k] = 0;
                } else if (filter_image[i][j][k] > _255) {
                    filter_image[i][j][k] = _255;
                }
            }
        }
    }
    // dezalocam imaginea
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            free(image[i][j]);
        }
        free(image[i]);
    }
    free(image);
    return filter_image;
}
