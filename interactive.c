#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "imageprocessing.h"
#include "bmp.h"
#define MAX 100
// aceasta functie schimba valorile stocate x si y;
void swap(int*x, int*y) {
    int temp = *x;
    *x = *y;
    *y = temp;
}
int main() {
    int ****Lista_imagini = NULL, *Lista_N = NULL, *Lista_M = NULL, *Filters_size = NULL;
    float ***Filetrs_list = NULL;
    int idx = 0, idxfil = 0;
    while (1) {
        char input[2];
        scanf("%s", input);
        if (!strcmp(input, "e")) {
            // daca input este egal cu "e" atunci trebuie sa dezalocam imaginile si filtrele care nu au fost sterse
            // la fel si vectorii in care sunt stocate marimile imaginilor si filtrelor
            for (int _i = 0; _i < idx; _i++) {
                for (int i = 0; i < Lista_N[_i]; i++) {
                    for (int j = 0; j < Lista_M[_i]; j++) {
                        free(Lista_imagini[_i][i][j]);
                    }
                    free(Lista_imagini[_i][i]);
                }
                free(Lista_imagini[_i]);
            }
            free(Lista_imagini);
            free(Lista_N);
            free(Lista_M);
            for (int i = 0; i < idxfil; i++) {
                for (int j = 0; j < Filters_size[i]; j++) {
                    free(Filetrs_list[i][j]);
                }
                free(Filetrs_list[i]);
            }
            free(Filetrs_list);
            free(Filters_size);
            // terminam executia ciclului si deoarece nu mai avem instructiuni dupa exit se termina si programul
            break;
        } else if (!strcmp(input, "l")) {
            // alocam si incarcam imaginea de dimensiune NxM aflată la calea path
            // pentru asta trebuie sa realocam vectorul de imagini si vectorii de marimi
            // marind dimensiunea lor cu un element
            int ****temporar_image = (int****)malloc((idx+1)*sizeof(int***));
            int*temporar_N = (int*)malloc((idx+1)*sizeof(int));
            int*temporar_M = (int*)malloc((idx+1)*sizeof(int));
            for (int i = 0; i < idx; i++) {
                temporar_image[i] = Lista_imagini[i];
                temporar_N[i] = Lista_N[i];
                temporar_M[i] = Lista_M[i];
            }
            free(Lista_imagini);
            free(Lista_N);
            free(Lista_M);
            Lista_imagini = temporar_image;
            Lista_N = temporar_N;
            Lista_M = temporar_M;
            int N = 0, M = 0;
            scanf("%d%d", &N, &M);
            char path[MAX];
            scanf("%s", path);
            Lista_imagini[idx] = (int***)malloc(N*sizeof(int**));
            for (int i = 0; i < N; i++) {
                Lista_imagini[idx][i] = (int**)malloc(M*sizeof(int*));
                for (int j = 0; j < M; j++) {
                    Lista_imagini[idx][i][j] = (int*)malloc(3*sizeof(int));
                }
            }
            Lista_N[idx] = N;
            Lista_M[idx] = M;
            read_from_bmp(Lista_imagini[idx], N, M, path);
            idx++;
        } else if (!strcmp(input, "s")) {
            // citim indexul la care vom salva imaginea
            int i = 0;
            scanf("%d", &i);
            char path[MAX];
            scanf("%s", path);
            // salvam imaginea de pe indexul respectiv
            write_to_bmp(Lista_imagini[i], Lista_N[i], Lista_M[i], path);
        } else if (!strcmp(input, "ah")) {
            // citim un index si aplicam functia flip_horizontal asupra imaginii de pe pozitia data
            int i = 0;
            scanf("%d", &i);
            Lista_imagini[i] = flip_horizontal(Lista_imagini[i], Lista_N[i], Lista_M[i]);
        } else if (!strcmp(input, "ar")) {
            // citim un index si aplicam functia rotate_left asupra imaginii de pe pozitia data
            int i = 0;
            scanf("%d", &i);
            Lista_imagini[i] = rotate_left(Lista_imagini[i], Lista_N[i], Lista_M[i]);
            // schimbam valorile la N si M dupa rotatie de pe pozitia respectiva
            // pentru a putea mai departe lucra cu imaginea modificata
            swap(&Lista_N[i], &Lista_M[i]);

        } else if (!strcmp(input, "ac")) {
            // citim un index si aplicam functia crop asupra imaginii de pe pozitia data
            int i = 0, x = 0, y = 0, h = 0, w = 0;
            scanf("%d%d%d%d%d", &i, &x, &y, &h, &w);
            Lista_imagini[i] = crop(Lista_imagini[i], Lista_N[i], Lista_M[i], x, y, w, h);
            // N si M trec in w si h dupa aplicarea functiei crop
            Lista_N[i] = w;
            Lista_M[i] = h;
        } else if (!strcmp(input, "ae")) {
            // citim un index si aplicam functia extend asupra imaginii de pe pozitia data
            int i = 0, rows = 0, cols = 0, new_R = 0, new_G = 0, new_B = 0;
            scanf("%d%d%d%d%d%d", &i, &rows, &cols, &new_R, &new_G, &new_B);
            Lista_imagini[i] = extend(Lista_imagini[i], Lista_N[i], Lista_M[i], rows, cols, new_R, new_G, new_B);
            // N si M trec in N + 2 * rows si M + 2 * cols dupa aplicarea functiei extend
            Lista_N[i]= 2 * rows + Lista_N[i];
            Lista_M[i]= 2 * cols + Lista_M[i];
        } else if (!strcmp(input, "ap")) {
            // aplicm operatia de paste cu parametrii dati imaginii de la indexul index_dst
            int d = 0, s = 0, x = 0, y = 0;
            scanf("%d%d%d%d", &d, &s, &x, &y);
            Lista_imagini[d] = paste(Lista_imagini[d], Lista_N[d], Lista_M[d], Lista_imagini[s],
                                   Lista_N[s], Lista_M[s], x, y);
         } else if (!strcmp(input, "cf")) {
            // initializam un filter si il alocam, dimensiunea lui fiind dim
            int dim = 0;
            scanf("%d", &dim);
            float***temp_filters = (float***)malloc((idxfil+1)*sizeof(float**));
            if (!temp_filters) {
                return -1;
            }
            for (int i = 0; i < idxfil; i++) {
                temp_filters[i] = Filetrs_list[i];
            }
            free(Filetrs_list);
            Filetrs_list = temp_filters;
            int*temp_filters_size = (int*)malloc((idxfil+1)*sizeof(int));
            if (!temp_filters_size) {
                return -1;
            }
            for (int i = 0; i < idxfil; i++) {
                temp_filters_size[i] = Filters_size[i];
            }
            free(Filters_size);
            Filters_size = temp_filters_size;
            Filters_size[idxfil] = dim;
            Filetrs_list[idxfil] = (float**)malloc(dim*sizeof(float*));
            if (!Filetrs_list[idxfil]) {
                free(Filetrs_list[idxfil]);
                return -1;
            }
            for (int i = 0; i < dim; i++) {
                Filetrs_list[idxfil][i] = (float*)malloc(dim*sizeof(float));
            }
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    scanf("%f", &Filetrs_list[idxfil][i][j]);
                }
            }
            idxfil++;
        } else if (!strcmp(input, "af")) {
            // aplicam filter-ul de la indexul ifil pe imaginea de la indexul iimg
            int iimg = 0, ifil = 0;
            scanf("%d%d", &iimg, &ifil);
            Lista_imagini[iimg] = apply_filter(Lista_imagini[iimg], Lista_N[iimg],
                                            Lista_M[iimg], Filetrs_list[ifil], Filters_size[ifil]);
        } else if (!strcmp(input, "df")) {
            // citim un index si stergem filter-ul de pe pozitia data
            int _i = 0;
            scanf("%d", &_i);
            for (int i = 0; i < Filters_size[_i]; i++) {
                free(Filetrs_list[_i][i]);
            }
            free(Filetrs_list[_i]);
            // stergem filter-ul de pe pozitia _i si mutam elementele care sunt dupa _i cu o pozitie in stanga
            for (int i = _i; i < idxfil - 1; i++) {
                Filetrs_list[i] = Filetrs_list[i+1];
                Filters_size[i] = Filters_size[i+1];
            }

            float***temp_fil = (float***)realloc(Filetrs_list, (idxfil-1)*sizeof(float**));
            Filetrs_list = temp_fil;
            int*temp_fil_size = (int*)realloc(Filters_size, (idxfil-1)*sizeof(int));
            Filters_size = temp_fil_size;
            idxfil--;
        } else if (!strcmp(input, "di")) {
            // citim un index si stergem imaginea de pe pozitia respectiva
            int _i = 0;
            scanf("%d", &_i);
            for (int i = 0; i < Lista_N[_i]; i++) {
                for (int j = 0; j < Lista_M[_i]; j++) {
                    free(Lista_imagini[_i][i][j]);
                }
                free(Lista_imagini[_i][i]);
            }
            free(Lista_imagini[_i]);
            // realocam vectorul Lista_imagini marimea fiind cu una mai mica decat curenta
            // dezalocam imaginea de pe indexul citit si mutam toate imaginile dupa cea stearsa cu o pozitie la stanga
            int****temporar_image = (int****)malloc((idx-1)*sizeof(int***));
            if (!temporar_image) {
                free(temporar_image);
                return -1;
            }
            // realocam vectorul Lista_N cu o marime mai mica cu 1 element
            // dezalocam marimea de pe pozitia _i si le mutam pe cele dupa _i cu o pozitie la stanga
            int*temporar_N = (int*)malloc((idx-1)*sizeof(int));
            if (!temporar_N) {
                free(temporar_N);
                return -1;
            }
            // realocam vectorul Lista_M cu o marime mai mica cu 1 element
            // dezalocam marimea de pe pozitia _i si le mutam pe cele dupa _i cu o pozitie la stanga
            int*temporar_M = (int*)malloc((idx-1)*sizeof(int));
            if (!temporar_M) {
                free(temporar_M);
                return -1;
            }
            for (int i = 0; i < _i; i++) {
                temporar_image[i] = Lista_imagini[i];
            }
            for (int i = _i; i < idx - 1; i++) {
                temporar_image[i] = Lista_imagini[i+1];
            }
            free(Lista_imagini);
            Lista_imagini = temporar_image;
            for (int i = 0; i < _i; i++) {
                temporar_N[i] = Lista_N[i];
            }
            for (int i = _i; i < idx - 1; i++) {
                temporar_N[i] = Lista_N[i+1];
            }
            free(Lista_N);
            Lista_N = temporar_N;
            for (int i = 0; i < _i; i++) {
                temporar_M[i] = Lista_M[i];
            }
            for (int i = _i; i < idx - 1; i++) {
                temporar_M[i] = Lista_M[i+1];
            }
            free(Lista_M);
            Lista_M = temporar_M;
            idx--;
        }
    }
    return 0;
}
