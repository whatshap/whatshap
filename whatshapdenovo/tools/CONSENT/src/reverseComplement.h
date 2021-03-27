#include <string>

#ifndef __REVCOMP__
#define __REVCOMP__

class rev_comp {
    public:
        static std::string run(std::string seq);

    protected:
        static rev_comp* _instance;
        
        static void build_instance();
        
        rev_comp();
            
        ~rev_comp();
        
    private:
        static char complement[int('t') + 1];
};

#endif