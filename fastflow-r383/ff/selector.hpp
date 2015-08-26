/* -*- Mode: C++; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */

/*! 

 *  \file stecilReduceOCL.hpp
 *  \ingroup high_level_patterns
 *
 *  \brief OpenCL map and non-iterative data-parallel patterns
 *
 */

/* ***************************************************************************
 *
 *  This program is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License version 3 as
 *  published by the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software Foundation,
 *  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 ****************************************************************************
 */

/*
 *  Authors:
 *    Massimo Torquati  (August 2015)
 *
 *  
 */

#ifndef FF_NODE_SELECTOR_HPP
#define FF_NODE_SELECTOR_HPP

namespace ff {


template<typename IN, typename OUT=IN>
class ff_nodeSelector: public ff_node_t<IN,OUT> {
protected:
    static bool ff_send_out_selector(void * task,
                                     unsigned long retry,
                                     unsigned long ticks, void *obj) {
        return reinterpret_cast<ff_node*>(obj)->ff_send_out(task, retry, ticks);
    }

public:

    ff_nodeSelector():selected(0) {}
    ff_nodeSelector(const IN &task):selected(0), inTask(const_cast<IN*>(&task)) {}
    
    // used to set tasks when running in a passive mode
    void setTask(const IN &task) { inTask = const_cast<IN*>(&task);  }

    void selectNode(size_t id) { selected = id; }

    int svc_init() { return nodeInit(); }

    OUT* svc(IN *in) {
        if (in == nullptr) {
            devices[selected]->svc(inTask);
            return ff_node_t<IN,OUT>::EOS;
        }
        OUT* out = (OUT*)(devices[selected]->svc(in));
        return out;
    }

    void svc_end() { nodeEnd(); }

    ff_node *getNode(size_t id) { 
        if (id >= devices.size()) return nullptr;
        return devices[id];
    }

    int nodeInit() {
        if (devices.size() == 0) return -1;
        for(size_t i=0;i<devices.size();++i)
            if (devices[i]->nodeInit()<0) return -1;
        return 0;
    }

    void nodeEnd() {
        for(size_t i=0;i<devices.size();++i)
            devices[i]->nodeEnd();
    }

    void addNode(ff_node &node) { 
        devices.push_back(&node); 
        node.registerCallback(ff_send_out_selector, this);                    
    }

    fftype getFFType() const   { 
        for(size_t i=0;i<devices.size();++i)
            if (devices[i]->getFFType() == OCL_WORKER) return OCL_WORKER;
        return WORKER;
    }

    int run(bool = false) { return ff_node::run();  }

    int wait() { return ff_node::wait(); }
    
    int run_and_wait_end() {
        if (run() < 0)	return -1;
        if (wait() < 0) return -1;
        return 0;
    }
            
protected:
    size_t      selected;
    IN         *inTask;
    std::vector<ff_node*> devices;    
};


} // namespace

#endif /* FF_NODE_SELECTOR_HPP */
