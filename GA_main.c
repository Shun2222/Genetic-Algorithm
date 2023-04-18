#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_GEN    50
#define POP_SIZE   128 //num of gene 2^n
#define LEN_CHROM  100  //length of gene
#define GEN_GAP    1.0
#define P_MUTATION 0.1
#define RANDOM_MAX 32767
#define BEFORE     0
#define AFTER      1
#define STARCOIN   20
#define SPACEMAX   42

int chrom[POP_SIZE][LEN_CHROM];
int fitness[POP_SIZE], maxfitness[MAX_GEN], avefitness[MAX_GEN];
int max, min, sumfitness;
int n_min, n_max;
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

void PrintFitnessHistory(){
    printf("--------------------Max fitness--------------------\n");
    printf("[");
    for(int i=0; i<MAX_GEN-1; i++){
        printf("%d, ", maxfitness[i]);
    }
    printf("%d]\n", maxfitness[MAX_GEN-1]);
    printf("--------------------Ave fitness--------------------\n");
    printf("[");
    for(int i=0; i<MAX_GEN-1; i++){
        printf("%d, ", avefitness[i]);
    }
    printf("%d]\n", avefitness[MAX_GEN-1]);
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
    n->nodeNum = -1;
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
        if(t->nodeNum == 32 && 43 <= num && num <= 53) {
            t = t->next2;
        } else {
            t = t->next1;
        }
    }
    return t;
}

void printlist(struct node *list) {
    int i;
    struct node *t;
    printf("{ ");
        for(i = 0; i <54; i++){
            t = list;
            printf("%d ", searchNode(t,i)->coin);
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

    return sugoroku->next1;
}

int ObjFunc(int i){
    int j, n, sumspace, star, sumcoin;
    struct node *space;
    
    star = sumspace = 0;
    sumcoin = 10;
    for(j=0; j<LEN_CHROM; j++){
        //さいころの値だけ次のマス
        space = searchNode(map, sumspace+chrom[i][j]);
        //スター獲得条件
        if(sumcoin>=STARCOIN && sumspace < 23 && sumspace+chrom[i][j] >= 23){
            star += 1; sumcoin -= STARCOIN;//スターの獲得とコインの消費
        } 
        sumcoin += space->coin;//コインの獲得
        sumspace += chrom[i][j];//進んだマスの更新
        if(sumspace==28) sumspace = 41;//強制移動マス
        if(sumspace>SPACEMAX){
            sumspace -= SPACEMAX;//一周
            sumcoin += 10;//一周完了報酬
        }
        sumcoin += 10;//毎ターン10コイン獲得
        
    }
    return star*1000  + sumcoin;//starが多いほど良い　starが同じ枚数ならcoinの枚数が多いほど良い　coinの枚数の最大でも3桁と考えた
}

void Statistics(int gen){
    int i;
    max = 0;
    min = fitness[0];
    sumfitness = 0;
    for(i=0; i<POP_SIZE; i++){
        if(fitness[i]>max){
            max = fitness[i];
            n_max = i;
        }
        if(fitness[i]<min){
            min = fitness[i];
            n_min = i;
        }
        sumfitness += fitness[i];
    }
    maxfitness[gen-1] = max;
    avefitness[gen-1] = sumfitness/POP_SIZE;
}

//Roulette Select 適応度に応じた割合
int RouletteSelect(){
    int i,sum, a=0;
    double rand;
    sum=0;
    rand = (double)Rand()/((double)(RANDOM_MAX+1));
    if(min<0) a=-min;
    for(i=0; i<POP_SIZE; i++){
        sum+=fitness[i]+a;
        if((double)sum/(double)(sumfitness+a*POP_SIZE)>rand) break;
    }
    if(i>127) i=127;
    return i;
}


int CheckSame(int *array, int len, int x){//arrayの0~len番目にxがあれば１なければ０ 
    if(len<=0) return 0;
    for(int i=0; i<len; i++){
        if(array[i]==x) return 1;
    }
    return 0;
}
int CheckMax(int *array, int start, int end){//arrayのstart~endまでにある最大
    int max = end;
    for(int i=start; i<end; i++){
        if(array[i]>array[max]) max = i;
    }
    return max;
}

void CopyArray(int *a, int *b){//aの配列にbの配列をcopy
    for(int i=0; i<POP_SIZE; i++){
        a[i] = b[i];
    }
}

int *Random_Pair(int *array, int n){
    int rand;
    int *result = malloc(sizeof(int) * n);
    int history[n];
    for(int i=0; i<n; i++){
        do{
            rand = Rand()%n;
        }while(CheckSame(history, i, rand));
        result[i] = array[rand];
        history[i] = rand;
    }
    //printf("\n");
    return result;
}

//Tournament Select 選択されたランダムな中から適応度が高い遺伝子
int TournamentSelect(int t){//nは2^i
    int i, rand, max=0;//tは一つのトーナメントに参加する遺伝子の数
    int participant[POP_SIZE];
    //トーナメント表作成
    for(i=0; i<t; i++){
        do{
            rand = Rand()%POP_SIZE;
        }while(CheckSame(participant, i, rand));//すでにparticipantにあったらやり直し
        participant[i] = rand;//participantに保存
        //printf("%d:%d ", rand, fitness[rand]); //参加者と適応度
        if(fitness[rand]>fitness[max]) max = rand;
    }
    //printf("\n max%d:%d\n", max, fitness[max]);//トーナメント勝者
    return max;
}

//あらかじめ選択されるものと数を決めるタイプ
//Ranking Select 適応度に依存しない
int Ranking(int n){//n番目に大きい適応度を返す
    int max=0;
    int copyfitness[POP_SIZE];
    CopyArray(copyfitness, fitness);
    for(int i=0; i<n; i++){
        max = CheckMax(copyfitness, 0, POP_SIZE-1);
        copyfitness[max] = 0;
    }
    return max;
}

int Calc_num(int x, int rank_max, int select_type){
    double a, b;
    int num, n=(int)(POP_SIZE*GEN_GAP);
    switch(select_type){
        case 3://比例
            a = -2.0*(double)n/(double)(rank_max*rank_max); b = 2.0*(double)n/(double)rank_max;
            num = (int)(a*((double)x-1.0)+b);
            break;
        case 4://反比例
            a = (double)(n)/log((double)(rank_max+1));
            num = (int)(a/(double)x);
            break;
        case 5://二次関数
            a = (double)(3*n)/(double)(rank_max*rank_max*rank_max), b=(double)rank_max;
            num = (int)(a*((double)x-b-1)*((double)x-b-1));
            break;
    }
    return num;
}

int *RankingSelect(int select_type){
    int i, j, k=0;
    int rank, rank_max=(int)((0.25)*(double)POP_SIZE);//rank:順位iの遺伝子のindex rank_max:選択される最低順位
    int num, n=(int)(POP_SIZE*GEN_GAP);//num:選択される回数　n:親の数
    if(n%2 != 0) n+=1;//ペアーを考慮
    if(n>POP_SIZE) n=POP_SIZE; //最大値の考慮
    int selected[n];//selectedが選択されたもの
    double a, b;

    for(i=1; (i<rank_max+1 || k<n); i++){
        rank = Ranking(i);//適応度の順位i番目
        num = Calc_num(i, rank_max, select_type);//選択個数　
        //printf("i:%d rank:%d n:%d k:%d\n", i, rank, n, k);
        //selectにrankをnum個保存
        for(j=0; j<num; j++){
            selected[k] = rank;
            k++;
            if(k==n) goto end;
        }
        //selectedが足りない場合の調整
        if(num<=0){
            selected[k] = rank;
            k++;
            if(k==n) goto end;     
        }
    }
end:
    return Random_Pair(selected, n);
}

//期待値選択
int *ExpectedValueSelect(){
    int i, j, k=0, a=0;
    int  n=(int)(POP_SIZE*GEN_GAP);
    if(n%2 != 0) n+=1;//ペアーを考慮

    int selected[n];
    double num, d, sum=0;
    if(min<0) a=-min;
    for(i=0; i<POP_SIZE; i++){
        d = ((double)(fitness[i]+a)/(double)(sumfitness+a*POP_SIZE));
        num = (double)d*n*5;
        //printf("num:%d d:%f fitness:%d\n", (int)num, d, fitness[i]);
        for(j=0; j<(int)num; j++){
            selected[k] = i;
            k++;
            if(k==n) goto next;
        }
    }
    for(i=1; n!=k; i++){
        selected[k] = Ranking(i);
        k++;
    }

next:
    return Random_Pair(selected, n);
}

//エリート選択 今回はGEN_GAP<1とすることで子供位に上書きさせないことで実現させる
void SaveElite(int n){
}





//交叉
void Crossover(int parent1, int parent2, int *child1, int *child2) {
    int min2;
    int n_cross;
    int i, j;

    /* 一番小さい値を子どもとしてセット */
    *child1 = n_min;
    /* 二番目に小さい値を見つける */
    min2 = max;
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
}

//突然変異
void Mutation(int child) {
    int n_mutate;
    double rand;

    rand = (double)Rand()/((double)(RANDOM_MAX+1));
    if(rand<P_MUTATION) {
        /* 突然変異位置 */
        n_mutate = Rand()%LEN_CHROM;
        /* 突然変異 */
        if(chrom[child][n_mutate] < 10) {  // 9以下なら+1する
            chrom[child][n_mutate] += 1;
        } else {                           // 10なら0にする
            chrom[child][n_mutate] = 1;
        }
        fitness[child] = ObjFunc(child);
    }
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
}


void Generation(int gen, int select){
    int parent1, parent2;
    int *parents;
    int child1, child2;
    int n_gen;
    int select_type = 1;
    int i;
    //データの表示
    Statistics(gen);
    //世代交代
    n_gen = (int)((double)POP_SIZE*GEN_GAP/2.0);

    if(select_type==0 || select_type==1){
        for(i=0; i<n_gen; i++){
            //Statistics(gen);
            //交叉する遺伝子の選択
            if(select_type==0){
                parent1 = RouletteSelect();
                parent2 = RouletteSelect();
            }else{
                parent1 = TournamentSelect(select);
                parent2 = TournamentSelect(select);
            }
            Crossover(parent1, parent2, &child1, &child2);//交叉
            Mutation(child1);   //child1 ある確率で突然変異 
            Mutation(child2);   //child2 ある確率で突然変異
        }
    }else{
        (select_type==2)? (parents=ExpectedValueSelect()) : (parents=RankingSelect(select_type));
        for(i=0; i<n_gen; i++){
            Statistics(gen);
        //交叉する遺伝子の選択
            parent1 = parents[2*i];
            parent2 = parents[2*i+1];
            Crossover(parent1, parent2, &child1, &child2);//交叉
            Mutation(child1);   //child1 ある確率で突然変異 
            Mutation(child2);   //child2 ある確率で突然変異
        }
    }
}

void test(int select_type){
    int gen;

    map = makeSugorokuMap();
    Initialize();

    for(gen=1; gen<=MAX_GEN; gen++) Generation(gen, select_type);
    printf("\n----------Results----------\n");
    printf("Max fitness:");
    PrintEachChromFitness(n_max);
    PrintFitnessHistory();
}

void main(int argc, char** argv){
    for(int i=1; i<11; i++){
        test(i);
    }
}
