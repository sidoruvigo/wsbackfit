interpret.sbformula <-
function(formula, method = "sback") {

    env <- environment(formula) 
    if(inherits(formula, "character"))          
        formula <- as.formula(formula)
    tf <- terms.formula(formula, specials = c("sb"))
    terms <- attr(tf, "term.labels")
    if(length(grep(":",terms))!=0)  stop("Symbol '*' is not allowed")
    
    nt <- length(terms)
   if(attr(tf, "response") > 0) {     
    	ns <- attr(tf, "specials")$sb - 1 # -1 for the response
    	response <- as.character(attr(tf, "variables")[2])
    } else {
    	ns <- attr(tf, "specials")$sb
    	response <- NULL
    }
    
    II <- list()
    h  <- list()
    partial = vector()
    k <- 0
    if(nt) {
        for (i in 1:nt) {
            if (i %in% ns) {
                k = k+1                   
                st <- eval(parse(text = terms[i]), envir = env)
                if(method ==  "sback" & st$cov[1] != "ONE") {
                	stop("The function sback does not allow for \"by\" variables")
                }
                II[[k]] <- st$cov
                h[[k]] <- st$h
                partial[k] <- terms[i]                                                          
            } else {
                k = k+1
                II[[k]]<- c("ONE", terms[i])
                h[[k]] <- 0
                partial[k] <- terms[i]
            }
        }           
    }       
    II <- if(length(II)) {
        matrix(unlist(II), nrow = 2)
    } else {
        matrix(0, nrow = 2)
    }       
    res <- list(response = response, II = II, h = unlist(h), npartial = k, partial = partial)
    res
}
