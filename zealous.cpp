#include<iostream>
#include<cmath>
#include<map>
#include<fstream>
#include<vector>
#include<string>
#include <iterator>
#include <random>
#include <algorithm>
#include <time.h>

using namespace std;

int main(int argc,char **argv){
    if(argc != 5){
        cerr << "You have to specify arguements to 'filename' 'epsilon' 'delta' 'm'" << endl;
        exit(1);
    }
    double epsilone = stod(argv[2]);
    double delta = stod(argv[3]);
    if(string(argv[4])!="o" and string(argv[4])!="m"){
        cerr << "The forth argument is need to be 'o' or 'm'." << endl;
        exit(1);
    }

    map<string,int>counter;

    const string filename=argv[1];
    ifstream reading_file;
    reading_file.open(filename,ios::in);
    if(reading_file.fail()){
        cerr<<"opening the file failed."<<endl;
        exit(-1);
    }

    while(!reading_file.eof()){
        string reading_line_buffer;
        getline(reading_file,reading_line_buffer);
        string user,domain;
        auto itr = reading_line_buffer.find(',');
        user=reading_line_buffer.substr(0,itr);

        counter[user]++;
    }
    reading_file.close();

    int num_of_user=0;
    long long total=0;
    for(auto itr = counter.begin(); itr != counter.end(); ++itr) {
        total+=itr->second;
        num_of_user++;
    }

    long long m;
    if(string(argv[4])=="o"){
        m=1;
    }else{
        m=total/num_of_user;
    }
    double lambda = 2*m/epsilone;

    int tau = 1;
    double taud,taud_next;
    do{
        taud = max(-lambda * log(2-2*exp(-1/lambda)),-lambda * log(2*delta/(num_of_user*m/tau)))+tau;
        tau++;
        taud_next = max(-lambda * log(2-2*exp(-1/lambda)),-lambda * log(2*delta/(num_of_user*m/tau)))+tau;
    }while(taud > taud_next);
    tau--;
    cerr << "m:" << m << " lambda:" << lambda << " tau:" << tau << " taud:" << taud << endl;


    ///////ファイル入力//////////
    map<string,vector<string>> users;

    reading_file.open(filename,ios::in);
    if(reading_file.fail()){
        cerr<<"opening the file failed."<<endl;
        exit(-1);
    }

    while(!reading_file.eof()){
        string reading_line_buffer;
        getline(reading_file,reading_line_buffer);
        string user,domain;
	    auto itr = reading_line_buffer.find(',');
        user=reading_line_buffer.substr(0,itr);
        domain = reading_line_buffer.substr(itr+1,reading_line_buffer.length()-itr-1);
        users[user].push_back(domain);
    }
    reading_file.close();

    ///////////各ユーザから最大m個のドメインをランダムに抽出してヒストグラムを構成////////////
    map<string,double> histogram;
    
    random_device seed_gen;// 乱数生成器を用意する
    default_random_engine engine {seed_gen()};

    for(map<string,vector<string>>::iterator itr = users.begin(); itr != users.end(); ++itr) {
        vector<string> selected_set,now;
        now = itr->second;
        sample(now.begin(),now.end(),back_inserter(selected_set),m,engine);
        
        for(int i=0;i<selected_set.size();i++){
            histogram[selected_set[i]]++;
        }
    }

    /////////ラプラスノイズ追加/////////////////
    exponential_distribution<>dist(lambda);
    vector<int> sign = {-1,1};
    
    vector<pair<double,string>>result;

    //////////閾値τ'以下の要素を削除//////////////
    for(map<string,double>::iterator itr = histogram.begin(); itr != histogram.end();++itr) {
        double freq = itr->second;
        vector<int> selected_sign;

        if(itr->second < tau)continue;
        sample(sign.begin(),sign.end(),back_inserter(selected_sign),1,engine);
        double noise = dist(engine)*selected_sign[0]/2;
        itr->second+=noise;
        
        if(itr->second >= max((double)tau,taud)){
            result.push_back(make_pair(itr->second,itr->first));
        }
    }

    sort(result.begin(),result.end());
    reverse(result.begin(),result.end());

    for(int i=0;i<result.size();i++){
        cout << result[i].first << " " << result[i].second << endl;
    }

    return 0;
}