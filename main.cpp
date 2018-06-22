#include <iostream>
#include <fstream>
#include <chrono>
#include <limits>

using namespace std;
int n_Nodes=100000;

class Node{ // Wezel rozwiazania i charakterystyczne dla niego wielkosci r[], q[], LB(K,c) (pomocnicze)
public:
    int* r;
    int* q;
    int priority;
    Node(int n): r(new int[n+1]),q(new int[n+1]) {}
    ~Node() {delete[] r; delete []q;}
};

void algorytmSchragePmtn(int &n, int *&r, int *&p, int *&q, long &cmax);
void algorytmSchrage(int &n, int *&o_perm, int *&r, int *&p, int *&q, long &cmax);

void heapSortDesc(int *&priority, int *&perm, int &n);
void heapSortAsc(int *&priority, int *&perm, int &n);


void heapPush(Node**& tabWez,int&n_heap);
Node* heapPop(Node**& tabWez,int&n_heap);

int Carlier(int &n, int*& perm, int *&r, int *&p, int *&q);

int main() {
    chrono::duration<double> czas_obliczen_NEH = chrono::duration<double>::zero();
    chrono::duration<double> elapsed_seconds;
    std::chrono::duration<double, std::nano> elapsed_microsec;
    int n;
    int *r,*p,*q, *perm;
    long cmax;

    string plik_wej = "carl.data.txt";
    string plik_wyj = "out_carl.txt";

    ifstream plik_in;
    ofstream plik_out;
    plik_out.open(plik_wyj, ofstream::trunc); plik_out.close(); // wyczysc plik wyjsciowy

    string linia;
    plik_in.open(plik_wej);

    while (!plik_in.eof()) {
        getline(plik_in, linia);
        if (linia[0] == 'd' && linia[1] == 'a') {
            plik_in >> n;

            perm = new int[n + 1];
            r = new int [n + 1];
            p = new int [n + 1];
            q = new int [n + 1];

            for (int i = 1; i <= n; ++i)
                plik_in >> r[i] >> p[i] >> q[i];

            chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
            cmax=Carlier(n,perm,r,p,q);
            chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
            elapsed_seconds = end - start;
            czas_obliczen_NEH += elapsed_seconds;

            plik_out.open(plik_wyj, ofstream::app);
            for (int i = 1; i <= n; ++i)
                plik_out << perm[i] << " ";
            plik_out << endl << cmax << endl;
            plik_out << "Czas obliczen: " << elapsed_seconds.count() << "s" << endl;
            plik_out.close();
            delete[] r; delete[] p; delete[] q; delete[] perm;
        }
    }
    plik_in.close();
    plik_out.open(plik_wyj, ofstream::app);
    plik_out << "Calkowity czas obliczen: " << czas_obliczen_NEH.count() << "s" << endl;
    plik_out.close();
//    cout << "Calkowity czas obliczen: " << czas_obliczen_NEH.count() << "s" << endl;

    return 0;
}

int Carlier(int &n, int*& o_perm, int *&r, int *&p, int *&q) { // TODO lista argumentow
    long u, LB, UB = std::numeric_limits<long>::max(); // dlugosc usz. Schrage, LowerBand, UpperBand
    int *C = new int [n+1]; C[0]=0; // terminy zakonczenia wykonywania zadan
    int a,b,c;  // indeksy zadan a, b, c
    int pi_i, pi_b, pi_c;   // zadania a, b, c
    int rk,pk,qk;   // R(K), P(K), Q(K)
    long hk, hk_c,LBA, LBB; // H(K), H(Kuc), LB(K,c) dla obu wezlow potomnych
    int *perm = new int[n+1];   // permutacja w danym wezle rozwiazania
    int temp;   // zmienna pomocnicza
    int n_heap=0;   // rozmiar kopca
    int ind;
    int elimin_threshold;
    Node* task; // biezacy wezel rozwiazania

    Node** taskHeap=new Node*[n_Nodes];  // tablica sluzaca jako kopiec rozwiazan oczekujacych
    for (int i=0; i<n_Nodes; ++i){
        taskHeap[i]=new Node(n+1);
    }
    for (int i=1; i<=n; ++i) {
        taskHeap[0]->r[i]=r[i];
        taskHeap[0]->q[i]=q[i];
    }
    task=taskHeap[0]; // wezel zerowy

    while( task != NULL) { // dopoki istnieje niezamkniety wezel

        c = 0; // c  domyslnie nieznalezione
        algorytmSchragePmtn(n, task->r, p, task->q, LB); // algorytm Schrage
        if (LB < UB) { // jesli zauktualizowany LB wiekszy lub rowny UB, zakoncz
            algorytmSchrage(n, perm, task->r, p, task->q, u); // mam cmax i permutacje

            for (int i = 1; i <= n; ++i) { // wyznacz wszystkie C, znajdz a, b
                pi_i = perm[i];
                if (task->r[pi_i] >= C[i - 1]) {  // wyznacz wszystkie C i znajdz a
                    a = i;
                    C[i] = task->r[pi_i] + p[pi_i];
                } else
                    C[i] = C[i - 1] + p[pi_i];

                if (C[i] + task->q[pi_i] == u) { // wyznacz b
                    b = i;
                    pi_b = perm[b];
                }
            } //mam a,b

            if (u < UB) { // aktualizacja UB, zapamietaj permutacje
                UB = u;
                for(int i=0; i<=n; ++i)
                    o_perm[i]=perm[i];
            }

            for (int i = b - 1, qb = q[pi_b]; i >= a; --i) // znajdz c
                if (task->q[perm[i]] < qb) {
                    c = i;
                    pi_c = perm[c];
                    break;
                }
        }
        // mam c
        if (c != 0) { // jesli nie znaleziono c (c=0), to zamknij wezel
            //znaleziono zadanie referencyjne
            // mam a,b,c, C dla (1,b)

            //wyznacz r(K),p(K),q(K),h(K)
            pk = C[b] - C[c];
            qk = task->q[pi_b];
            rk = std::numeric_limits<int>::max();
            for (int i = b; i > c; --i)
                rk = min(task->r[perm[i]], rk);

            hk = rk + pk + qk;

            // wezel 1A (zadanie c poprzedza zbior K ) - q'
            temp = max(task->q[pi_c], pk + qk);//korekta q, zadanie c przed blokiem
            hk_c = pk + p[pi_c] + min(qk, temp) + min(rk, task->r[pi_c]);
            LB = max(LB, hk);
            LBA = max(LB, hk_c); // oblicz  dolne ograniczenie dla nastepnika A
            if (LBA < UB) { // dodaj zadanie, jesli dolne ograniczenie mniejsze od UB
                // pomocnicze dolne ograniczenie jako priorytet w kopcu
                taskHeap[n_heap + 1]->priority = LBA;
                for (int i = 1; i <= n; ++i) {//przepisz r i q
                    taskHeap[n_heap + 1]->q[i] = task->q[i];
                    taskHeap[n_heap + 1]->r[i] = task->r[i];
                }
                taskHeap[n_heap + 1]->q[pi_c] = temp; // aktualizacji q'
                //dodatkowe testy eliminacyjne
                elimin_threshold=UB-hk;
                for (int i = 1; i <= n; ++i) // dla zbioru J\K
                    if (i < c || b < i) {
                        ind = perm[i];
                        if (p[ind] > elimin_threshold) {
                            if (task->r[ind] + p[ind] + pk + qk >= UB)
                                taskHeap[n_heap + 1]->r[ind] = max(task->r[ind], rk + pk);
                            if (rk + p[ind] + pk + task->q[ind] >= UB)
                                taskHeap[n_heap + 1]->q[ind] = max(task->q[ind], qk + pk);
                        }
                    }

                heapPush(taskHeap, n_heap); // dodanie zadania do kopca
            }

            // wezel 1A (zadanie c na koncu zbioru K ) - r'
            temp = max(task->r[pi_c], rk + pk); //korekta r, zadanie c za blokiem
            hk_c = pk + p[pi_c] + min(qk, task->q[pi_c]) + min(rk, temp);
            LBB = max(LB, hk_c); // oblicz  dolne ograniczenie dla nastepnika B
            if (LBB < UB) { // dodaj zadanie, jesli dolne ograniczenie mniejsze od UB
                // pomocnicze dolne ograniczenie jako priorytet w kopcu
                taskHeap[n_heap + 1]->priority = LBB;
                for (int i = 1; i <= n; ++i) {//przepisz r i q
                    taskHeap[n_heap + 1]->q[i] = task->q[i];
                    taskHeap[n_heap + 1]->r[i] = task->r[i];
                }
                taskHeap[n_heap + 1]->r[pi_c] = temp; // akutalizacja r'
                //dodatkowe testy eliminacyjne
                elimin_threshold=UB-hk;
                for (int i = 1; i <= n; ++i) // dla zbioru J\K
                    if (i < c || b < i) {
                        ind = perm[i];
                        if (p[ind] > elimin_threshold) {
                            if (task->r[ind] + p[ind] + pk + qk >= UB)
                                taskHeap[n_heap + 1]->r[ind] = max(task->r[ind], rk + pk);
                            if (rk + p[ind] + pk + task->q[ind] >= UB)
                                taskHeap[n_heap + 1]->q[ind] = max(task->q[ind], qk + pk);
                        }
                    }
                heapPush(taskHeap, n_heap); // dodanie zadania do kopca
            }
        }
        // pobierz kolejne zadanie z kopca, jesli istnieje i jego priorytet jest niewiekszy od UB
        do{
            task = heapPop(taskHeap, n_heap);
        }while(task != NULL && task->priority >= UB);
    }
    //zwolnienie pamieci
    delete[] C;
    for (int i=0; i<n_Nodes; ++i){
        delete taskHeap[i];
    }
    delete []taskHeap;
    delete []perm;

    return UB; //zwroc dlugosc uszeregowania
}


void algorytmSchragePmtn(int &n, int *&r, int *&p, int *&q, long &cmax) {
    int t = 0, e;
    int *N = new int[n + 1], *G = new int[n + 1];
    int *p_pom = new int[n + 1];
    int nG = 0, nN = n;
    int l=0;
    cmax=0;
    q[0]= std::numeric_limits<int>::max();;

    for (int j = 1; j <= n; ++j){
        N[j] = j;
        p_pom[j]=p[j];
    }
    heapSortDesc(r, N, nN);

    while (nG || nN) {
        while (nN && (r[N[nN]] <= t)) { // jesli element o najmniejszym r
            e = N[nN];
            G[++nG] = e;
            --nN; // zmniejsz rozmiar i automatycznie wyklucz najmniejszy element
            if(q[e]>q[l]){
                p_pom[l]=t-r[e];
                t=r[e];
                if (p_pom[l]>0)
                    G[++nG]=l;
            }
        }
        if (nG == 0)
            t = r[N[nN]];
        else {
            heapSortAsc(q, G, nG);
            e = G[nG--];
            l=e;
            t += p_pom[e];
            cmax = max(cmax, (long) t + q[e]);
        }
    }
    delete[] N; delete[] G; delete[] p_pom;
}
void algorytmSchrage(int &n, int *&perm, int *&r, int *&p, int *&q, long &cmax) {
    int t = 0, k = 0, e;
    int *N = new int[n + 1], *G = new int[n + 1];
    int nG = 0, nN = n;
    cmax=0;

    for (int j = 1; j <= n; ++j) N[j] = j;
    heapSortDesc(r, N, nN);

    while (nG || nN) {
        while (nN && (r[N[nN]] <= t)) { // jesli element o najmniejszym r
            e = N[nN];
            ++nG;
            G[nG] = e;
            --nN; // zmniejsz rozmiar i automatycznie wyklucz najmniejszy element
        }
        if (nG == 0)
            t = r[N[nN]];
        else {
            heapSortAsc(q, G, nG);
            e = G[nG--];
            perm[++k] = e;
            t += p[e];
            cmax = max(cmax, (long) t + q[e]);
        }
    }
    delete[] N; delete[] G;
}

void heapSortDesc(int *&priority, int *&perm, int &n) {
    int p_ind;  // rodzic (parent index)
    int ch_ind; //  dziecko (child index)
    int gch_ind; // wieksze dziecko (greater child index)
    int x, y;

    for (int i = 2; i <= n; i++) { // budowanie kopca
        ch_ind = i;
        p_ind = ch_ind / 2;
        x = priority[perm[i]];
        y = perm[i];

        while ((p_ind > 0) && ((priority[perm[p_ind]] > x)|| ((x == priority[perm[p_ind]]) &&(perm[p_ind]>y)))) {
            perm[ch_ind] = perm[p_ind];
            ch_ind = p_ind;
            p_ind = (ch_ind >> 1);
        }
        perm[ch_ind] = y;
    }


    for (int i = n; i > 1; i--)      // rozebranie kopca
    {  // element najwiekszy na koniec
        swap(perm[1], perm[i]);
        p_ind = 1;
        ch_ind = 2;
        while (ch_ind < i) //z
        {
            if ((ch_ind + 1 < i) &&
                ((priority[perm[ch_ind + 1]] < priority[perm[ch_ind]]) ||
                 ((priority[perm[ch_ind + 1]] == priority[perm[ch_ind]]) && (perm[ch_ind+1] < perm[ch_ind])))) //jesli istnieje prawe dziecko i jest mniejsze od lewego
                gch_ind = ch_ind + 1; // mniejsze prawe
            else
                gch_ind = ch_ind;   // mniejsze lewe
            if ((priority[perm[gch_ind]] > priority[perm[p_ind]]) || ((priority[perm[gch_ind]] == priority[perm[p_ind]]) &&(perm[p_ind]<perm[gch_ind])))
                break; // jesli mniejsze dziecko jest wieksze od rodzica przerwij 'while'
            swap(perm[p_ind], perm[gch_ind]);// jesli nie, zamien rodzica z wiekszym dzieckiem
            p_ind = gch_ind;
            ch_ind = (p_ind << 1);
        }
    } // zmniejsz rozmiar kopca
}
void heapSortAsc(int *&priority, int *&perm, int &n) {
    int p_ind;  // rodzic (parent index)
    int ch_ind; //  dziecko (child index)
    int gch_ind; // wieksze dziecko (greater child index)
    int x, y;

    for (int i = 2; i <= n; i++) { // budowanie kopca
        ch_ind = i;
        p_ind = ch_ind / 2;
        x = priority[perm[i]];
        y = perm[i];

        while ((p_ind > 0) && ((priority[perm[p_ind]] < x)|| ((x == priority[perm[p_ind]]) &&(perm[p_ind]<y)))) {
            perm[ch_ind] = perm[p_ind];
            ch_ind = p_ind;
            p_ind = (ch_ind >> 1);
        }
        perm[ch_ind] = y;
    }


    for (int i = n; i > 1; i--)      // rozebranie kopca
    {  // element najwiekszy na koniec
        swap(perm[1], perm[i]);
        p_ind = 1;
        ch_ind = 2;
        while (ch_ind < i) //z
        {
            if ((ch_ind + 1 < i) &&
                ((priority[perm[ch_ind + 1]] > priority[perm[ch_ind]]) ||
                 ((priority[perm[ch_ind + 1]] == priority[perm[ch_ind]]) && (perm[ch_ind+1] > perm[ch_ind])))) //jesli istnieje prawe dziecko i jest mniejsze od lewego
                gch_ind = ch_ind + 1; // mniejsze prawe
            else
                gch_ind = ch_ind;   // mniejsze lewe
            if ((priority[perm[gch_ind]] < priority[perm[p_ind]]) || ((priority[perm[gch_ind]] == priority[perm[p_ind]]) &&(perm[p_ind]>perm[gch_ind])))
                break; // jesli mniejsze dziecko jest wieksze od rodzica przerwij 'while'
            swap(perm[p_ind], perm[gch_ind]);// jesli nie, zamien rodzica z wiekszym dzieckiem
            p_ind = gch_ind;
            ch_ind = (p_ind << 1);
        }
    } // zmniejsz rozmiar kopca
}


void heapPush(Node**& tabWez,int&n_heap) {
    int ch,p;
    n_heap++;
    if(n_heap>=n_Nodes) {
        cerr << "Przerost kopca" << endl;
        exit(-2);
    }

    Node* temp = tabWez[n_heap];
    int prior= temp->priority;

    ch = n_heap;
    p = ch / 2;

    while(p > 0 && tabWez[p]->priority > prior) {
        tabWez[ch] = tabWez[p];
        ch = p;
        p = ch / 2;
    }

    tabWez[ch] = temp;
}
Node *heapPop(Node **&tabWez, int &n_heap) {
    if (!n_heap)
        return NULL;

    int p, ch;
    Node *tmp;
    Node *tmp2; // latka

    tmp = tabWez[n_heap];
    tmp2 = tabWez[0];
    tabWez[0] = tabWez[1];
    tabWez[n_heap] = tmp2;
    --n_heap;
    if(n_heap) {
        p = 1;
        ch = 2;

        while (ch <= n_heap) {
            if (ch + 1 <= n_heap && tabWez[ch + 1]->priority < tabWez[ch]->priority) ch++;
            if (tmp->priority <= tabWez[ch]->priority) break;
            tabWez[p] = tabWez[ch];
            p = ch;
            ch = 2 * ch;
        }
        tabWez[p] = tmp;
    }

    return tabWez[0];
}

