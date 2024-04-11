# 
function getm(s::Any)
        """
        grabs methods
        asin|>getm  
        ?asin
        @code_llvm readf|>getm|>first 
        """
        methods(s);
    end

    