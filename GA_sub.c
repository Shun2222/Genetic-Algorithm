#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_GEN    30
#define POP_SIZE   128 //num of gene 2^n
#define LEN_CHROM  10  //length of gene
#define GEN_GAP    0.2
#define P_MUTAION  0.1
#define RANDOM_MAX 32767
#define BEFORE     0
#define AFTER      1
#define STARCOIN   20
#define ELITENUM   5


int chrom[POP_SIZE][LEN_CHROM];
int fitness[POP_SIZE], elite[ELITENUM];
int max, min, sumfitness;
int n_min;
struct node *map;

static unsigned long int next = 1;

int Rand(void){
    next = next*1103515245+12345;
    return (unsigned int)(next/65536)%32768;
}

void Srand(unsigned int seed){
    next = seed;
}

//i番目の遺伝子データの表示
void PrintEachChromFitness(int i){
    int j;
    printf("[%d]", i);
    for(j=0; j<LEN_CHROM; j++){
        printf("%d ", chrom[i][j]);  //i番目の遺伝子の表示
    }
    printf(":%d\n", fitness[i]);    //適合度の表示
}

void PrintChromFitness(){
    int i;
    for(i=0; i<POP_SIZE; i++) PrintEachChromFitness(i); //すべての遺伝子データと適合度の表示
}

//学習経過の表示
void PrintStatistics(int gen){
    printf("[gen=%2d] max=%d min=%d sumfitness=%d ave=%f\n", gen, max, min, sumfitness, (double)sumfitness/(double)POP_SIZE);
}


void PrintCrossover(int flag, int parent1, int parent2, int child1, int child2, int n_cross){
    switch(flag){
        case BEFORE:
            printf("parent1 |");PrintEachChromFitness(parent1);
            printf("parent2 |");PrintEachChromFitness(parent2);
            printf("delete1 |");PrintEachChromFitness(child1);
            printf("delete2 |");PrintEachChromFitness(child2);
            printf("n_cross=%d\n", n_cross);
            break;
        case AFTER:
            printf("child1 |");PrintEachChromFitness(child1);
            printf("child2 |");PrintEachChromFitness(child2);
            printf("--------------------\n");
            break;
    }
}


void PrintMutation(int flag, int child, int n_mutate){
    switch(flag){
        case BEFORE:
            printf("child(OLD)|"); PrintEachChromFitness(child);
            printf("n_mutate=%d\n", n_mutate);
            break;
        case AFTER:
            printf("child(NEW)|"); PrintEachChromFitness(child);
            printf("--------------------\n");
            break;
    }
}

struct node {
    int star;
    int coin;
    int turn;
    int nodeNum;
    struct node *next1;
    struct node *next2;
};

struct node *initlist() {
    struct node *n;
    n = (struct node *) malloc (sizeof (struct node));
    n->next1 = NULL;
    n->next2 = NULL;
    return n;
}

void nodeAppendNext1(struct node *p, int num) {
    struct node *n, *tmp;
    n = (struct node*) malloc (sizeof(struct node));
    tmp = p;
    while(tmp->next1 != NULL) {
        tmp = tmp->next1;
    }
    n->nodeNum = num;
    n->next1 = NULL;
    n->next2 = NULL;
    tmp->next1 = n;
}

void nodeAppendNext2(struct node *p, int num) {
    struct node *n, *tmp;
    n = (struct node*) malloc (sizeof(struct node));
    tmp = p;
    while(tmp->next2 != NULL) {
        tmp = tmp->next2;
    }
    n->nodeNum = num;
    n->next2 = NULL;
    tmp->next2 = n;
}

struct node *searchNode(struct node *list, int num) {
    struct node *t;
    t = list;
    while(t->nodeNum != num) {
        t = t->next1;
    }
    return t;
}

void printlist(struct node *p)
{
    int i;
    printf("{ ");
        for(i = 0; i < 100; i++){
        printf("%d ", p->nodeNum);
        if(p->nodeNum == 32 || p->nodeNum == 53) {
            p = p->next2;
        }
        else {
        p = p->next1;
        }
        }
    
   printf("}\n");
}



void setCoin(struct node *list) {
    int coin[54] = {
        0,10,3,0,0,3,-3,0,0,3,10,3,0,0,3,
        3,0,0,3,10,10,3,3,3,3,3,10,3,0,3,
        3,-3,10,3,10,10,3,0,0,-3,0,3,-20,
        -3,-3,-3,3,10,3,3,10,10,0,-20
    };

    struct node *t;
    for(int i = 0; i < 54; i++) {
        t = list;
        t = searchNode(t,i);
        t->coin = coin[i];
    }
}


struct node *makeSugorokuMap() {
    struct node *sugoroku;
    sugoroku = initlist();

    for(int i = 0; i <= 42; i++) {
        nodeAppendNext1(sugoroku, i);
    }

    struct node *t1, *t2;
    t1 = searchNode(sugoroku, 1); t2 = searchNode(sugoroku, 42);
    t2 -> next1 = t1;

    t1 = searchNode(sugoroku, 19); t2 = searchNode(sugoroku, 32);
    nodeAppendNext2(searchNode(sugoroku, 32), 43);
    t2 = t2->next2;
    for(int i = 44; i <= 53; i++) {
        nodeAppendNext1(t2, i);
    }
    t2 = searchNode(t2, 53);
    t2->next1 = t1;

    t1 = searchNode(sugoroku, 28); t2 = searchNode(sugoroku, 41);
    t1->next2 = t2;

   setCoin(sugoroku);

    return sugoroku;
}

//基本の目的関数
/*
int ObjFunc(int i){
    int j;
    int count=0;
    for(j=0; j<LEN_CHROM; j++){
        if(chrom[i][j] == 1) count++;
    }
    return count;
}

void CopyArray(int *a, int *b){
    for(int i=0; i<POP_SIZE; i++){
        for(int j=0; j<LEN_CHROM; j++){
            a[i][j] = b[i][j];
        }
    }
}

Survival(int gen){
    int survivor;
    int noob[LEN_CHROM], nextgen[POP_SIZE][LEN_CHROM];
    int n_gen;
    int i;
    //データの表示
    Statistics();
    PrintStatistics(gen);

    //世代交代
    n_gen = (int)(POP_SIZE/2);
    for(i=0; i<n_gen; i++){
        survivor = RouletteSelect();
        nextgen[i] = chrom[survivor];
    }
    for(i=n_gen; i<POP_SIZE; i++){
        for(j=0; j<LEN_CHROM; j++){
            nextgen[i][j] = Rand()%10+1; 
        }
    }
}
*/

int ObjFunc(int i){
    int j, n, sumspace, star, sumcoin;
    struct node *space = map;
    
    star = sumspace = sumcoin = 0;
    for(j=0; j<LEN_CHROM; j++){
        //さいころの値だけ次のマス
        for(n=0; n<chrom[i][j]; n++){
            space = space->next1;//次のマス
        }
        sumcoin += space->coin;//コインの獲得
        //スター獲得条件
        if(sumcoin>STARCOIN && sumspace < 23 && sumspace+chrom[i][j] > 23){
            star += 1; sumcoin -= STARCOIN;//スターの獲得とコインの消費
        } 
        sumspace += chrom[i][j];//進んだマスの更新
    }
    return star + sumcoin* 0.001;//starが多いほど良い　starが同じ枚数ならcoinの枚数が多いほど良い　coinの枚数の最大でも3桁と考えた
}

void Statistics(){
    int i;
    max = 0;
    min = POP_SIZE;
    sumfitness = 0;
    for(i=0; i<POP_SIZE; i++){
        if(fitness[i]>max) max = fitness[i];
        if(fitness[i]<min){
            min = fitness[i];
            n_min = i;
        }
        sumfitness += fitness[i];
    }
}

//Roulette Select 適応度に応じた割合
int RouletteSelect(){
    int i,sum;
    double rand;
    sum=0;
    rand = (double)Rand()/((double)(RANDOM_MAX+1));
    for(i=0; i<POP_SIZE; i++){
        sum+=fitness[i];
        if((double)sum/(double)sumfitness>rand) break;
    }
    return i;
}
/*//まとめて
int *RouletteSelect(int n){
    int i;
    int sum;
    double rand;
    int *results = malloc(sizeof(int) * n);
    for(i=0; i<n; i++){
        sum=0;
        rand = (double)Rand()/((douoble)(RANDOM_MAX+1));
        for(j=0; j<POP_SIZE; j++){
            sum+=fitness[j];
            if((double)sum/(double)sumfitness>rand) break;
        }
        results[i] = j;
    }
    return results;
}

//Tournament Select 選択されたランダムな中から適応度が高い遺伝子
bool CheckSame(int *array, int len, int x){
    if(len<=0) return false;
    for(int i=0; i<len; i++){
        if(array[i]==x) return true;
    }
    return false;
}
int CheckMax(int *array, int start, int end){
    int max = end;
    for(int i=start; i<end; i++){
        if(array[i]>array[max]) max = i;
    }
    return max;
}
//ひとつ
int TournamentSelect(){//nは2^i
    int i, j, x, rand, max=0, t=10;//tは一つのトーナメントに参加する遺伝子の数
    int history[POP_SIZE];
    //トーナメント表作成
    for(j=0; j<t; j++){
        do{
            rand = Rand()%LEN_CHROM;
        }while(CheckSame(history, j, rand));
        history[j] = rand;
        if(fitness[rand]>fitness[max]) max = rand;
    }
    return max;
}
//まとめて
int *TournamentSelect(int n){//nは2^i
    int i, j, x, rand, max=0, t=10;//tは一つのトーナメントに参加する遺伝子の数
    int history[POP_SIZE];
    int *results = malloc(sizeof(int) * n);
    //トーナメント表作成
    for(i=0; i<n; i++){
        max = 0;
        for(j=0; j<t; j++){
            do{
                rand = Rand()%LEN_CHROM;
            }while(CheckSame(history, j, rand));
            history[j] = rand;
            if(fitness[rand]>fitness[max]) max = rand;
        }
        results[i] = max;
    }
    return results;
}

//あらかじめ選択されるものと数を決めるタイプ
//Ranking Select 適応度に依存しない
int *RankingSelect(int n){
    int m = {}; 
    int i;
    int *x = malloc(sizeof(int) * n);
    int *results = malloc(sizeof(int) * n);
    for(i=0; i<n; i++){
         results[i] = CheckMax(fittness, 0, POP_SIZE-1);
         x[i] = fitness[results[i]]; fitness[results[i]] = 0;
    }
    for(i=0; i<n; i++){
        fitness[results[i]] = x[i];
    }
    return results;
}

void SaveElite(int n){
    for(i=0; i<ELITENUM; i++){
         elite[i] = CheckMax(fittness, 0, POP_SIZE-1);
         fitness[elite[i]] = 0;
    }
}

//期待値選択
int ExpectedValueSelect(int n){
    int i,sum;
    double rand;
    sum=0;
    rand = (double)Rand()/((douoble)(RANDOM_MAX+1));
    for(i=0; j<POP_SIZE; i++){
        sum+=fitness[i];
        if((double)sum/(double)sumfitness>rand) break;
    }
    return i;
}
*/






//交叉
void Crossover(int parent1, int parent2, int *child1, int *child2) {
    int min2;
    int n_cross;
    int i, j;

    /* 一番小さい値を子どもとしてセット */
    *child1 = n_min;
    /* 二番目に小さい値を見つける */
    min2 = POP_SIZE;
    for(i = 0; i < POP_SIZE; i++) {
        if(i != *child1) {
            if(min <= fitness[i] && fitness[i] < min2) {
                min2 = fitness[i];
                *child2 = i;
            }
        }
    }

    /* 交叉位置 */
    n_cross = Rand()%(LEN_CHROM-1)+1;  /* n_cross = 1, ..., 9 */

    /* 交叉 */
    PrintCrossover(BEFORE, parent1, parent2, *child1, *child2, n_cross);
    for(j = 0; j < n_cross; j++) {
        chrom[*child1][j] = chrom[parent1][j];
        chrom[*child2][j] = chrom[parent2][j];
    }
    for(j = n_cross ; j < LEN_CHROM; j++) {
        chrom[*child1][j] = chrom[parent2][j];
        chrom[*child2][j] = chrom[parent1][j];
    }
    fitness[*child1] = ObjFunc(*child1);
    fitness[*child2] = ObjFunc(*child2);
    PrintCrossover(AFTER, parent1, parent2, *child1, *child2, n_cross);
}

//突然変異
void Mutation(int child) {
    
}


//初期値はランダム
void Initialize(){
    int i,j;
    for(i=0; i<POP_SIZE; i++){
        for(j=0; j<LEN_CHROM; j++){
            chrom[i][j] = Rand()%10+1; 
        }
        fitness[i] = ObjFunc(i);
    }
    printf("First Population\n");
    PrintChromFitness();
    printf("--------------------\n");
}

void Generation(int gen){
    int parent1, parent2;
    int child1, child2;
    int n_gen;
    int i;
    //データの表示
    Statistics();
    PrintStatistics(gen);

    //世代交代
    n_gen = (int)((double)POP_SIZE*GEN_GAP/2.0);
    for(i=0; i<n_gen; i++){
        Statistics();
        //交叉する遺伝子の選択
        parent1 = RouletteSelect();
        parent2 = RouletteSelect();
        Crossover(parent1, parent2, &child1, &child2);//交叉
        //Mutation(child1);   //child1 ある確率で突然変異 
        //Mutation(child2);   //child2 ある確率で突然変異
    }
}

void main(int argc, char** argv){
    int gen;
    map = makeSugorokuMap();
    Initialize();
    for(gen=1; gen<=MAX_GEN; gen++) Generation(gen);
}


/* Memo
遺伝子の長さ = ターン数 10
!マス 爆弾設置　コイン払う　or 強制移動 
VS かけたコインget 
クッパマーク　ランダム　
クローバー　ルーレット　アイテムorコイン get
分岐　ルーレット4つ　外れリセット　あたりそのまま
star coin 
*/

