;;;(setq Info-directory "/yale/lib/emacs/info")
(setq auto-mode-alist
      (quote (("\\.c$" . c++-mode)
	      ("\\.h$" . c++-mode)
	      ("\\.cc$" . c++-mode)
	      ("\\.hh$" . c++-mode) 
	      ("\\.cpp$" . c++-mode)
	      ("\\.hpp$" . c++-mode)
              ("\\.java$" . c++-mode)
	      ("\\.f$" . fortran-mode)
	      ("\\.tex$" . LaTeX-mode)
	      ("\\.lisp$" . lisp-mode)
	      ("\\.el$" . emacs-lisp-mode)
	      ("\\.Z$" . uncompress-while-visiting)
	      ("\\.F$" . unfreeze-while-visiting)
	      ("\\.z$" . gunzip-while-visiting)
	      ("\\.t$" . t-mode)
	      ("\\.m$" . matlab-mode)
	      ("\\.txt$" . text-mode))))
(display-time)
(put 'upcase-region 'disabled nil)
(set-background-color "black")
(set-foreground-color "white")
(set-cursor-color "white")
(global-set-key "" (quote goto-line))
(global-set-key [(meta c)] 'kill-ring-save)
(display-time)
(setq c-brace-offset -2)
(setq c-label-offset 0)
(autoload (quote ispell-buffer) "ispell" "Interactive speller." t)
(autoload (quote uncompress-while-visiting) "uncompress" nil t)
(if (eq window-system (quote x)) (global-set-key "" (quote backward-delete-char)))
(setq visible-bell t)
(setq inhibit-startup-message t)
(autoload (quote uncompress-while-visiting) "uncompress" nil t)
(autoload (quote unfreeze-while-visiting) "unfreeze" nil t)
(autoload (quote gunzip-while-visiting) "gunzip" nil t)
(autoload 'c++-mode "c++-mode"   nil t)

;; transparency
(set-frame-parameter nil 'alpha 80) 
(setq default-frame-alist
      (append
       (list
	'(active-alpha . 0.9)  ;; active frame
	'(inactive-alpha . 0.8) ;; non active frame
	) default-frame-alist)
      )
(tool-bar-mode)
;(mac-key-mode)

(defun ifdef-comment ()
  "Insert an ifdef comment at point."
  (interactive)
  (insert "#ifdef documentation\n")
  (insert "=========================================================================\n")
  (insert "\n")
  (insert "     program: ")
  (insert (buffer-name))
  (insert "\n          by: justin gardner\n")
  (insert "        date: ")
  (shell-command "date +%D" t)
  (exchange-point-and-mark)
  (insert "\n")
  (insert "=========================================================================\n")
  (insert "#endif\n")
  (set-mark (point)))

(defun matlab-stub ()
  "Insert matlab function stub."
  (interactive)
  (insert "% ")
  (setq commandname (substring (buffer-name) 0 (- (length (buffer-name)) 2)))
  (insert (buffer-name))
  (insert "\n%\n")
  (insert "%        $Id:$ \n")
  (insert "%      usage: ")
  (insert commandname)
  (insert "()\n%         by: justin gardner\n")
  (insert "%       date: ")
  (shell-command "date +%D" t)
  (exchange-point-and-mark)  
  (insert "%    purpose: \n")
  (insert "%\n")
  (insert "function retval = ")
  (insert commandname)
  (insert "()\n")
  (insert "\n")
  (insert "% check arguments\nif ~any(nargin == [1])\n")
  (insert "  help ")
  (insert commandname)
  (insert "\n  return\nend\n\n")
  (set-mark (point)))

(defun box-comment (string)
  "Boxed comment for matlab programs."
  (interactive "s")
  (let ((len (+ (length string) 10))
	(i 0))
    (while (< i len)
      (insert "%")
      (setq i (+ i 1)))
    (insert "\n")
    (insert "%    ")
    (insert string)
    (insert "    %\n")
    (setq i 0)
    (while (< i len)
      (insert "%")
      (setq i (+ i 1))))
  (insert "\n"))

(defun line-comment (string)
  "Boxed comment for matlab programs."
  (interactive "s")
  (insert "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  (insert "%    ")
  (insert string)
  (insert "\n")
  (insert "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;   by: justin gardner
;; date: 12/6/2002
;;
;; these functions are used to generate sequential numbers
;; they use two global variables, count which contains
;; the current number and count_0 which is set to the 
;; number of digits you want printed (i.e. if count_0 
;; equals 2 than you will get 00,01,02...). You can
;; restart counting by typing M-x count_start_at. To start
;; at a specific number, then type C-u num M-x count_start_at
;; Each time you want to insert a number into the buffer
;; run M-x count_next
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defvar count 1)
(defvar count_0 2)
(defun count_start ()
  "starts the count for count_next"
  (interactive)
  (setq count 1))
(defun count_start_at (num)
  "starts the count at num for count_next"
  (interactive "p")
  (setq count num))
(defun count_zeros (num)
  "Sets the number of zeros to print out for count_next"
  (interactive "p")
  (setq count_0 num))
(defun count_next ()
  "Inserts the next number into the buffer.
   Use (Ctrl-u num M-x count_start) to set first number
   or use (M-x count_start_at) to start at 1
   Use count_zeros to set number of digits to print"
  (interactive)
  (insert (format (format "%%0%dd" count_0) count))
  (setq count (+ count 1))
  (set-mark (point)))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(setq auto-mode-alist
      (quote (("\\.c$" . c++-mode)
	      ("\\.h$" . c++-mode)
	      ("\\.cc$" . c++-mode)
	      ("\\.hh$" . c++-mode) 
	      ("\\.cpp$" . c++-mode)
	      ("\\.hpp$" . c++-mode)
              ("\\.java$" . c++-mode)
	      ("\\.f$" . fortran-mode)
	      ("\\.tex$" . LaTeX-mode)
	      ("\\.lisp$" . lisp-mode)
	      ("\\.el$" . emacs-lisp-mode)
	      ("\\.Z$" . uncompress-while-visiting)
	      ("\\.F$" . unfreeze-while-visiting)
	      ("\\.z$" . gunzip-while-visiting)
	      ("\\.t$" . t-mode)
	      ("\\.m$" . matlab-mode)
	      ("\\.txt$" . text-mode))))
	     

(setq matlab-mode-hook 'turn-off-auto-fill)
(setq c-mode-hook 'turn-off-auto-fill)
;(setq LaTeX-mode-hook 'turn-on-auto-fill)
;(setq emacs-lisp--mode-hook 'turn-on-auto-fill)
;(setq lisp-mode-hook 'turn-on-auto-fill)
;(setq c++-mode-hook 'turn-on-auto-fill)
;(setq c-mode-hook 'turn-on-auto-fill)
;(setq default-major-mode 'text-mode)
(defun my-matlab-mode-hook ()
  (setq fill-column 500))		; where auto-fill should wrap

(put 'eval-expression 'disabled nil)


;
; Tempus nachempfunden:
;
(defvar fliesstext nil "Merker fuer fliesstext-mode")
(make-variable-buffer-local 'fliesstext)
(setq minor-mode-alist (cons '(fliesstext " Fliess") minor-mode-alist))
(setq auto-mode-alist (cons '("^/tmp/Re" . fliesstext-mode) auto-mode-alist))
(setq auto-mode-alist (cons '("^/tmp/xmail" . fliesstext-mode) auto-mode-alist))

(defun turn-on-fliesstext-mode ()
  "Schaltet den Fliesstextmodus ein."
  (interactive "*")
  (setq fliesstext nil)
  (fliesstext-mode))

(defun fliesstext-mode () 
  "Schaltet zwischen normalem und Fliesstextmodus um. 
In letzterem reformatiert ein Leerzeichen automatisch den laufenden Absatz."
  (interactive "*")
  (setq fliesstext (not fliesstext))
  (if fliesstext
      (progn
	(local-set-key " " 'my-space)
	(local-set-key "\M-q" 'my-fill-paragraph))
    (progn
      (local-set-key " " 'self-insert-command)
      (local-set-key "\M-q" 'fill-paragraph))))

(defun my-space ()
  "Reformatiert den laufenden Absatz und fuegt ein Leerzeichen ein.
Wird verwendet im fliesstext-mode."
  (interactive "*")
  (my-fill-paragraph nil)
  (insert " "))

(defun my-fill-paragraph (arg)
  "Bricht Paragraph um.
Mit Praefixargument auch Randausgleich."
  (interactive "P")
  (save-excursion
    (beginning-of-line)
    (if (not (looking-at paragraph-separate))
	(progn
	  (if (not (looking-at paragraph-start))
	      (beg-of-par nil))
	  (if (looking-at paragraph-separate)
	      (next-line 1))
	  (let ((end (point)))
	    (end-of-par nil)
	    (fill-region-as-paragraph (point) end arg))))))

(defun beg-of-par (arg)
  (interactive "P")
  (if (re-search-backward paragraph-start nil t arg)
      (goto-char (match-beginning 0))
    (goto-char (point-min))))

(defun end-of-par (arg)
  (interactive "P")
  (if (looking-at paragraph-start)
      (forward-char 1))
  (if (re-search-forward paragraph-start nil t arg)
      (goto-char (match-beginning 0))
    (goto-char (point-max))))


;; -*- Mode: Emacs-Lisp -*-
;; File:            c++-mode.el
;; Description:     Mode for editing C++ code
;; Authors:         1992 Barry A. Warsaw, Century Computing Inc.
;;                  1987 Dave Detlefs  (dld@cs.cmu.edu)
;;                   and Stewart Clamen (clamen@cs.cmu.edu)
;;                  Done by fairly faithful modification of:
;;                  c-mode.el, Copyright (C) 1985 Richard M. Stallman.
;; Last Modified:   $Date: 1992/05/08 20:47:38 $
;; Version:         $Revision: 2.40 $

(defun comment-region ()
  "Comment out all lines in a region between mark and current point by
inserting comment-start in front of each line."
  (interactive)
  (let* ((m      (if (eq (mark) nil) (error "Mark is not set!") (mark)))
	 (start  (if (< (point) m) (point) m))
	 (end    (if (> (point) m) (point) m))
	 (mymark (copy-marker end)))
    (save-excursion
	(goto-char start)
	(while (< (point) (marker-position mymark))
	    (beginning-of-line)
	    (insert comment-start)
	    (beginning-of-line)
	    (forward-line 1)))))

(defun uncomment-region ()
  "Uncomment all lines in region between mark and current point by deleting
the leading comment-start from each line, if any."
  (interactive)
  (let* ((m      (if (eq (mark) nil) (error "Mark is not set!") (mark)))
	 (start  (if (< (point) m) (point) m))
	 (end    (if (> (point) m) (point) m))
	 (mymark (copy-marker end))
	 (len    (length comment-start))
	 (char   (string-to-char comment-start)))
    (save-excursion
	(goto-char start)
	(while (< (point) (marker-position mymark))
	    (beginning-of-line)
	    (if (looking-at (concat " *" comment-start))
		(progn
		  (zap-to-char 1 char)
		  (delete-char len)))
	    (beginning-of-line)
	    (forward-line 1)))))

(set-fill-column 80)

(cond ((fboundp 'global-font-lock-mode)
       ;; Turn on font-lock in all modes that support it
       (global-font-lock-mode t)
       ;; Maximum colors
       (setq font-lock-maximum-decoration t)))

(cond ((fboundp 'global-font-lock-mode)
       ;; Customize face attributes
       (setq font-lock-face-attributes
             ;; Symbol-for-Face Foreground Background Bold Italic Underline
             '((font-lock-comment-face       "cyan")
               (font-lock-string-face        "red")
               (font-lock-keyword-face       "green")
               (font-lock-function-name-face "red")
               (font-lock-variable-name-face "purple")
               (font-lock-type-face          "yellow")
               ))
       ;; Load the font-lock package.
       (require 'font-lock)
       ;; Maximum colors
       (setq font-lock-maximum-decoration t)
       ;; Turn on font-lock in all modes that support it
       (global-font-lock-mode t)))

(custom-set-variables
  ;; custom-set-variables was added by Custom.
  ;; If you edit it by hand, you could mess it up, so be careful.
  ;; Your init file should contain only one such instance.
  ;; If there is more than one, they won't work right.
 '(custom-file nil))
(custom-set-faces
  ;; custom-set-faces was added by Custom.
  ;; If you edit it by hand, you could mess it up, so be careful.
  ;; Your init file should contain only one such instance.
  ;; If there is more than one, they won't work right.
 )
(load "~/proj/matlab/matlab.el")
(setq default-frame-alist '(
                             (border-color          . "green")
                             (foreground-color      . "white")
                             (background-color      . "black")
                             (vertical-scroll-bars  . nil)
                             (cursor-color          . "yellow")
                             (cursor-type           . box)
                             (top . 1) (left . 270) (width . 120) (height . 90)))

(message "Hello, friend")
