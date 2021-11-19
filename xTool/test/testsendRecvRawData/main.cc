/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
// xTool includes
#include "xMPIEnv.h"
#include "xRawDataExchanger.h"
//system include
#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <algorithm>
#include <functional>
#include <random>
#include <string>
#include <fstream>
#include <chrono>

int main(int argc, char *argv[]){
    xtool::xMPIEnv::init(argc,argv);
    int max_message_size = 10;
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    //select number of proc to send to
    //unsigned seed =  std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(111123);
    std::uniform_int_distribution<int> nb_target_distribution(0,mpi_size);
    int nbtargets = nb_target_distribution(generator);

    //select targets
    std::set<int> targets;
    std::uniform_int_distribution<int> target_distribution(0,mpi_size-1);
    while (targets.size() != nbtargets){
          targets.insert(target_distribution(generator));
    }
    std::uniform_int_distribution<int> message_size_distribution(0,max_message_size);
    std::uniform_real_distribution<double> message_elem_distribution(0.,1.);

    std::map<int, std::vector<double> > messages;
    //targets.clear();
    if (mpi_rank == 0) targets.insert(1);
    for(auto to : targets){
        std::vector<double> message_to( message_size_distribution(generator));
        std::generate(message_to.begin(), message_to.end(),
                      [&generator, &message_elem_distribution]()
                      {return message_elem_distribution(generator);});
        messages[to] = std::move(message_to);
    }
    for  (auto to_message :messages ){
        int to = to_message.first;
        const auto & message = to_message.second;
        std::ofstream file("proc_"+std::to_string(mpi_rank)+"_send_to_"+std::to_string(to));
        file << message.size() << std::endl;
        for (auto d : message) file << d << " ";
        file << std::endl;
    }
    std::function<void (int, const std::vector<double> & )> recvfunc
            = [&mpi_rank](int from, const std::vector<double> &message ){
        std::ofstream file("proc_"+std::to_string(mpi_rank)
                           +"_receive_from_"+std::to_string(from));
        file << message.size() << std::endl;
        for (auto d : message) file << d << " ";
        file << std::endl;
    };

    xtool::sendRecvRawData(messages, recvfunc, MPI_COMM_WORLD);

    return xtool::xMPIEnv::finalize();
}

